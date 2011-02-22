
## Copyright (c) 2001-2009, Scott D. Peckham
## January 2009  (converted from IDL)
## August 2009
## October 2009  (routines to allow more output file formats)

#-------------------------------------------------------------------
#  Functions:
#
#      write_profile()
#      Number_of_Samples()
#      save_step()
#      As_Grid()
#      Pixel_Var()
#      Stack()
#      Profile_Var()

#      add_file_extension()      # (12/8/09)
#      open_new_gs_file()
#      convert_scalar_to_grid()  # (8/18/09)
#      save_as_grid_to_file()    # (8/18/09)
#          OR "save_to_grid_stack()" OR "add_to_grid_stack()"  
#      close_gs_file()

#      open_new_ts_file()
#      write_ts_file_header()
#      write_ts_file_line()
#      close_ts_file()

#      write_pixel_file_header()   # (OBSOLETE)
#      write_pixel_file_line()     # (OBSOLETE)

#-------------------------------------------------------------------
from numpy import *
import numpy

import sys  # (for sys.byteorder)

from tf_utils import *

import rts_files
import rti_files
import ncgs_files
# import bov_files
import rtg_files

#-------------------------------------------------------------------
def write_profile(file_unit, ptr, outlet_IDs, nz, tmstr):

    import idl_func
    
    #--------------------------------------------------------
    # If (*ptr) is a 1D profile, and therefore the same
    # for all outlet_IDs, then the Profile_Var function
    # now returns a 1D array.  Otherwise, it returns a
    # 2D array, with a separate profile for each outlet_ID.
    #--------------------------------------------------------
    var = Profile_Var(ptr, outlet_IDs)
    nd  = numpy.ndim(var)

    
    n   = size(outlet_IDs)  # (needed below)  #####
    
    #----------------
    # For debugging
    #----------------
    #print 'nz = ', nz
    #print 'n  = ', n
    #print 'size(var) = ', size(var)
    
    file_unit.write(tmstr)
    
    if (nd == 1):    
        format = '%18.10f'
        for j in xrange(nz):
            file_unit.write( format % var[j] )
    else:    
        f = '(' + str(n) + 'F18.10)'
        for j in xrange(nz):
            out_str = idl_func.string(var[j,:], format=f)
            file_unit.write( out_str )
    
#   write_profile()
#-------------------------------------------------------------------
def Number_of_Samples(stop_vars, sample_dt, dt):

    #-----------------------------------------------
    # T_Stop_Model_Estimate uses max basin area to
    # get a dynamic upper bound for T_stop_model.
    # T_stop_model = (t_Pstop + (L_max / v_avg))
    #-----------------------------------------------
    # durations    = *pv.durations
    # Qp_fraction  = stop_vars.Qp_fraction
    # basin_areas  = (*grid_vars.basin_areas)[0]
    # t_stop_model = T_Stop_Model_Estimate(basin_areas, durations, $
    #                                     Qp_fraction)
    
    #------------------------------
    # Get number of sampled times
    #------------------------------
    I2PY_expr = stop_vars.method
    
    if (I2PY_expr == 0):    
        #-----------------------------------
        # Run until Q drops to P% of Qpeak
        #-----------------------------------
        t_stop_model = stop_vars.t_stop_model
        n_samps = int32(t_stop_model * float64(60) / sample_dt) + int32(1)
        
        #***********************************
        #4/19/06.  EXPERIMENTAL
        #***********************************
        #** n_samps = (n_samps * 4L)
        
        #*** print,'t_stop_model = ', t_stop_model
        #*** print,'n_samps      = ', n_samps
        #*** n_samps = 10000L   ;*******
        
        #secs_per_year = 31449600L   ;(60 * 60 * 24 * 365)
        #n_samps = long(secs_per_year / sample_dt) + 1L
    elif (I2PY_expr == 1):    
        #-----------------------------------
        # Run for a specified time (model)
        #-----------------------------------
        t_stop_model = stop_vars.t_stop_model
        n_samps = int32(t_stop_model * float64(60) / sample_dt) + int32(1)
    elif (I2PY_expr == 2):    
        #--------------------------------------
        # Run for a specified number of steps
        #--------------------------------------
        n_steps = stop_vars.n_steps
        t_stop_model = n_steps * dt
        n_samps = int32(t_stop_model / sample_dt) + int32(1)
        #--------------------------------------
        # Could also update t_stop_model when
        # n_steps is first set by user
        #--------------------------------------
        #** t_stop_model = stop_vars.t_stop_model
        #** n_samps = long(t_stop_model * 60d / sample_dt) + 1L
    else:
        raise RuntimeError('no match found for expression')
    
    return n_samps
    
#   Number_of_Samples()
#-------------------------------------------------------------------
def save_step(sample_dt, main_dt):

    step = ceil((maximum(sample_dt, main_dt)) / chan_dt).astype('Int32')
    
    return step
    
#   save_step()
#-------------------------------------------------------------------
def As_Grid(var, nx, ny):

    if (size(var) > 1):    
        return float32(var)
    
    #---------------------------------
    # Var is scalar; convert to grid
    #---------------------------------
    grid = zeros([ny, nx], dtype='Float32')
    return (grid + float32(var))
    
#   As_Grid()
#-------------------------------------------------------------------
def Pixel_Var(var, outlet_IDs):
    
    #---------------------------------
    # Is variable a grid or scalar ?
    #---------------------------------
    ### if (numpy.size(var) > 1):
    if (numpy.rank(var) > 0):
        return float32(var[outlet_IDs])
        ## return float32(var.flat[outlet_IDs])
    else:
        #-------------------------------------------------------
        # (3/16/07) Bug fix.  This gets used in case of q0,
        # which is a scalar when INFIL_ALL_SCALARS is true.
        # Without this, don't get a value for every outlet_ID.
        #-------------------------------------------------------
        # NB!  outlet_rows = outlet_IDs[0] as long as
        #      outlet_IDs is a numpy.where-style tuple.
        #-------------------------------------------------------                
        n_outlets = numpy.size(outlet_IDs[0])
        vector = zeros(n_outlets, dtype='Float32')
        return (vector + float32(var)) 
    
#   Pixel_Var()
#-------------------------------------------------------------------
def Stack(var, nx, ny):

    #---------------------
    # Get the dimensions
    #---------------------
    ndim = numpy.ndim(var)
    
    if (size(ndim) == 3):    
        #-------------------------
        # Variable is a 3D array
        #-------------------------
        return float32(var)
    else:    
        #-----------------------------------------
        # Var is 1D profile; convert to 3D array
        #-----------------------------------------
        nz = size(var)
        stack = zeros([nz, ny, nx], dtype='Float32')
        for k in xrange(nz):
            stack[k,:,:] = float32(var[k])
        return stack
    
#   Stack()
#-------------------------------------------------------------------
def Profile_Var(var, outlet_IDs):

    #---------------------
    # Get the dimensions
    #---------------------
    ndims = numpy.ndim(var)
    
    if (ndims == 1):    
        #---------------------------------
        # Variable is a 1D profile, and
        # is the same for all outlet_IDs
        #---------------------------------
        return float32(var)
    else:    
        #---------------------------------
        # Variable is a 3D array; return
        # a profile for each outlet_ID
        #---------------------------------
        nz    = size(var,0)  # (RE-CHECK IF Z-INDEX IS 0.)
        n_IDs = size(outlet_IDs)
        profiles = zeros([nz, n_IDs], dtype='Float32')
        for k in xrange(nz):
            layer = var[k]
            vals = layer[outlet_IDs]
            profiles[k,:] = float32(vals)
        return profiles
    
#   Profile_Var()
#-------------------------------------------------------------------
def add_file_extension( file_name, extension='.nc' ):

    n_ext = len(extension)
    if (file_name[-n_ext:].lower() != extension.lower()):
        return (file_name + extension)
    else:
        return file_name

#   add_file_extension()
#-------------------------------------------------------------------
def open_new_gs_file(file_name, info=None,
                     ## format='RTG',    # (switched off, 12/02/09)
                     format='NCGS',
                     dtype='float32',
                     grid_name='X',
                     long_name='Unknown',
                     units_name='None',
                     MAKE_BOV=True,
                     nx=None, ny=None, dx=None, dy=None):

    #---------------------------
    # Was grid info provided ?
    #---------------------------
    if (info != None):
        info.file_name = file_name
        info.data_type = rti_files.get_rti_data_type( dtype )
        if (nx != None): info.ncols = nx
        if (ny != None): info.nrows = ny
        if (dx != None): info.xres  = dx
        if (dy != None): info.yres  = dy
    else:
        if (nx != None) and (ny != None) and \
           (dx != None) and (dy != None):
            info = rti_files.make_info( file_name, nx, ny, dx, dy )
        else:
            print 'ERROR during open_new_gs_file().'
            print '      Grid info not provided.'
            print ' '
            return -1
        
    #--------------------------------
    # Try to get info from RTI file
    #--------------------------------
##    RTI_file = rti_files.try_to_find_rti_file( file_name )
##    if (RTI_file != 'none'):
##        info = rti_files.read_info( RTI_file )
##        info.file_name = file_name
##        info.data_type = rti_files.get_rti_data_type( dtype )
##        if (nx != None): info.ncols = nx
##        if (ny != None): info.nrows = ny
##        if (dx != None): info.xres  = dx
##        if (dy != None): info.yres  = dy
##    else:
##        if (nx != None) and (ny != None) and \
##           (dx != None) and (dy != None):
##            #------------------------
##            # Create a new RTI file
##            #------------------------
##            info = rti_files.make_info( file_name, nx, ny, dx, dy )
##            rti_files.write_info( file_name, info )
##        else:
##            print 'ERROR: Grid info not provided and cannot'
##            print '       find RTI file in working directory.'
##            print ' '
##            return

    #-----------------------------------------
    # Open new file of specified file format
    #-----------------------------------------
    if (format == 'RTS') :
        rts = rts_files.rts_file()
        rts_file_name = add_file_extension(file_name, '.rts')
        OK  = rts.open_new_file( rts_file_name, info, grid_name,
                                MAKE_BOV=MAKE_BOV )
        if (OK): return rts
    #----------------------------------------------------------
    elif (format == 'NCGS'):
        ncgs = ncgs_files.ncgs_file()
        nc_file_name = add_file_extension(file_name, '.nc')
        OK = ncgs.open_new_file( nc_file_name, info,
                                 grid_name, long_name,
                                 units_name )
        if (OK): return ncgs
    #----------------------------------------------------------
    elif (format == 'RTG') :
        rtg = rtg_files.rtg_file()
        rtg_file_name = add_file_extension(file_name, '.rtg')
        OK = rtg.open_new_file( rtg_file_name, info, grid_name,
                                ADD_INDEX=True,
                                MAKE_BOV=MAKE_BOV )
        if (OK):
            rtg.file_name_base = file_name  ######
            return rtg
    #----------------------------------------------------------
    else:
        OK = False
        print 'SORRY, "' + format + '" output is not supported.'
        print ' '

    if not(OK):
        print 'SORRY, Could not open new grid stack file:'
        print '   ' + file_name
        return -1
    
#   open_new_gs_file()
#-------------------------------------------------------------------
def convert_scalar_to_grid(scalar, nx, ny):

    grid = zeros([ny, nx], dtype='Float32')
    
    return (grid + float32(scalar))

#   convert_scalar_to_grid()
#-------------------------------------------------------------------
def save_as_grid_to_file(file_unit, var, var_name, nx, ny):
    
    #-----------------------------------------------------
    # Note that this function currently converts data to
    # 4-byte float, and should probably have a "dtype"
    # argument to be more general.
    #-----------------------------------------------------
    MAKE_BOV = True
    format = file_unit.format

    #-------------------------------------------------------
    # Don't do next line, because it changes byte order of
    # var such that write_pixel_file_line() has problems.
    #-------------------------------------------------------    
    ## if (SWAP_ENDIAN): var.byteswap(True)
    var  = float32(var)  ###########
    var2 = var.copy()
    var2 = float32(var2)   ######### (Float64 to Float32)

    #---------------------------------------------------------
    # Byteswap if necessary to maintain consistent byteorder
    # within a set of RTG and RTS files.  This is not needed
    # for netCDF (NCGS and NCG below).
    #---------------------------------------------------------
    if (file_unit.info.SWAP_ENDIAN) and (format in ['RTG','RTS']):
        var2.byteswap(True)
    
    if (numpy.rank(var2) == 0):
        grid = convert_scalar_to_grid(var2, nx, ny)  #(returns float32)
    else:
        grid = var2   # (reference synonym)

    #-----------------------------------------
    # Write grid to file in specified format
    #-----------------------------------------
    if  (format == 'RTS'):
        file_unit.add_grid( grid )
    elif (format == 'NCGS'):
        file_unit.add_grid( grid, var_name )
    elif (format == 'RTG'):
        file_unit.write_grid( grid )  # (includes close())
        OK = file_unit.open_new_file( file_unit.file_name_base,
                                      info=None,  # (Need this here ?)
                                      var_name=var_name,
                                      MAKE_RTI=False,
                                      MAKE_BOV=MAKE_BOV,
                                      ADD_INDEX=True) ######
        if not(OK): print 'ERROR while opening new RTG file.'
    elif (format == 'NCG'):
        file_unit.add_grid( grid, var_name )
        #--------------------------------------------
        # Close file with 1 grid and open a new one
        #--------------------------------------------
        file_unit.close()
        suffix   = str(file_unit.time_index).zfill(5)
        new_file = file_unit.file_name + suffix
        nx = file_unit.nx
        ny = file_unit.ny
        dx = file_unit.dx
        dy = file_unit.dy
        long_name  = file_unit.long_name
        units_name = file_unit.units_name
        OK = file_unit.open_new_file( new_file, nx, ny, dx, dy,
                                      var_name, long_name,
                                      units_name )
        if not(OK): print 'ERROR while opening new NCG file.'
    else:
        pass

#   save_as_grid_to_file()
#-------------------------------------------------------------------
def close_gs_file(file_unit):

    file_unit.close_file()

#   close_gs_file()
#-------------------------------------------------------------------
#-------------------------------------------------------------------
def open_new_ts_file(ts_file_name, var_name, IDs,
                     time_units='minutes'):

    try:
        ts_unit = open(ts_file_name, 'w')
        write_ts_file_header(ts_unit, IDs, var_name=var_name,
                         time_units=time_units)
        return ts_unit
    except:
        print 'SORRY, Could not open new time series file:'
        print '   ' + ts_file_name
        return -1
    
#   open_new_ts_file()
#-------------------------------------------------------------------
def write_ts_file_header(file_unit, outlet_IDs, var_name='F',
                         time_units='minutes', SILENT=True):

    #------------------------------------------------------
    # Notes:  This is currently hardwired to have columns
    #         that are each 15 characters wide.
    #------------------------------------------------------
    rows = outlet_IDs[0]
    cols = outlet_IDs[1]
    if not(SILENT):
        print 'In model_output.write_ts_file_header():'
        print 'outlet_rows =', rows
        print 'outlet_cols =', cols
        print ' '
        
    #---------------------------------------------------------
    # NB!  Allow "outlet_IDs" to be passed as either a tuple
    #      (numpy.where style indices) or as a 1D array
    #      (calendar style indices).
    #---------------------------------------------------------
##    if ('tuple' in str(type(outlet_IDs))):
##        rows = outlet_IDs[0]
##        cols = outlet_IDs[1]
##        print 'In model_output.write_ts_file_header():'
##        print 'outlet_rows =', rows
##        print 'outlet_cols =', cols
##        print ' '
##    else:
##        # Requires that we pass in nx
##        rows = (outlet_IDs / nx)
##        cols = (outlet_IDs % nx)
        
    n_outlets = size(rows)
    col_width = 15

    ustr_map = {'minutes':'min', 'seconds':'sec', \
                'hours':'hrs', 'days':'days', \
                'months':'mon', 'years':'yrs'}
    time_str = 'Time [' + ustr_map[time_units.lower()] + ']'
    
    #--------------------------
    # Create the heading text
    #--------------------------
    file_unit.write( time_str.rjust(col_width) )
    for k in xrange(n_outlets):
        str1 = str(rows[k])
        str2 = str(cols[k])
        str3 = var_name + '[' + str1 + ',' + str2 + ']'
        file_unit.write( str3.rjust(col_width) ) 
    file_unit.write("\n")
    
    #------------------------------------
    # Draw a horizontal line of hyphens
    #------------------------------------
    width = (n_outlets + 1) * col_width
    hline = ''.ljust(width, '-')
    file_unit.write(hline + "\n")

#   write_ts_file_header()
#-------------------------------------------------------------------
def write_ts_file_line(file_unit, time, var, outlet_IDs):

    #---------------------------------------------------------
    # NB!  Allow "outlet_IDs" to be passed as either a tuple
    #      (numpy.where style indices) or as a 1D array
    #      (calendar style indices).

    #      Units of time are not specified on purpose, but
    #      are available from a component's state.
    #---------------------------------------------------------

    #--------------------------------------------------------
    # Assume that outlet_IDs are tuples as with numpy.where
    # Note: ('tuple' in str(type(outlet_IDs))) is True.
    #--------------------------------------------------------
    rows = outlet_IDs[0]
    cols = outlet_IDs[1]

    #------------------------------------------
    # Assume that outlet_IDs is a 1D array of
    # long-integer "calendar-style" indices
    #------------------------------------------
##        rows = (outlet_IDs / nx)
##        cols = (outlet_IDs % nx)
        
    n_outlets = size(rows)
    col_width = 15

    #--------------------------------------------------------
    # If "var" is a grid, subscript with outlet_ID to get a 1D
    # array of values.  If "var" is scalar, return a vector
    # with the scalar value repeated once for each "outlet".
    #--------------------------------------------------------    
    vals = Pixel_Var(var, outlet_IDs)
    
    tstr = ('%15.7f' % time)
    file_unit.write( tstr.rjust(col_width)  )
    for k in xrange(n_outlets):
        vstr = ('%15.7f' % vals[k])   ## (same as v. 1.5b)
        file_unit.write( vstr.rjust(col_width) ) 
    file_unit.write("\n")

#   write_ts_file_line()
#-------------------------------------------------------------------
def close_ts_file(file_unit):
    
    file_unit.close_file()

#   close_ts_file()
#-------------------------------------------------------------------
#-------------------------------------------------------------------
##def save_as_grid_to_file_OLD(file_unit, var, nx, ny, SWAP_ENDIAN,
##                         long_name, units_name,
##                         file_type='NC'):
##                         #####   file_type='RTS'):
##
##    #------------------------------------------------------
##    # Don't do this, because it changes the byte order of
##    # var such that write_pixel_file_line() has problems.
##    #------------------------------------------------------
##    ## if (SWAP_ENDIAN): var.byteswap(True)
##    var  = float32(var)  ###########
##    var2 = var.copy()
##    var2 = float32(var2)   ######### (Float64 to Float32)
##    
##    if (SWAP_ENDIAN): var2.byteswap(True)
##    
##    if (numpy.rank(var2) == 0):
##        grid = convert_scalar_to_grid(var2, nx, ny)  #(returns float32)
##    else:
##        grid = var2  # (reference synonym)
##
##    if (file_type == 'RTS'):
##        grid.tofile(file_unit)
##    elif (file_type == 'NC'):
##        write_grid_to_netcdf(file_unit, grid, var_name,
##                             long_name, units_name)
##
###   save_as_grid_to_file_OLD()
#-------------------------------------------------------------------
##def open_new_netcdf_file(nc_file, nx, ny,
##                         var_name, long_name, units_name):
##
##    try:
##        import Nio
##        print 'Imported Nio version: ' + Nio.__version__
##    except:
##        python_version = sys.version[:3]
##        print 'SORRY, Cannot write netCDF files because'
##        print 'the "Nio" package cannot be imported.'
##        print 'Note that "PyNIO" is only installed for'
##        print 'Python version 2.6 on "beach".'
##        print 'The current Python version is:' + python_version
##        print ' '
##        return
##
##    #----------------------------
##    # Does file already exist ?
##    #--------------------------------------------
##    # Overwrite not allowed; must remove first.
##    #--------------------------------------------
##    if (os.path.exists( nc_file )):
##        print 'Deleting existing file: ' + nc_file
##        os.remove( nc_file )
##
##    #----------------------------
##    # Does file already exist ?
##    #----------------------------
####    if (os.path.exists(nc_file)):
####        print 'WARNING: Found a file named:'
####        print '  ' + nc_file
####        print 'in the current working directory.'
####        print 'Are you sure you want to overwrite it?'
####        ## os.remove(nc_file)   # (delete existing file)
####        return
##        
##    #-------------------------------------
##    # Open a new netCDF file for writing
##    #-------------------------------------
##    # Sample output from time.asctime():
##    #     "Thu Oct  8 17:10:18 2009"
##    #-------------------------------------
##    opt = Nio.options()
##    opt.PreFill = False            # (for efficiency)
##    opt.HeaderReserveSpace = 4000  # (4000 bytes, for efficiency)
##    history = "Created by TopoFlow 3.0 on "
##    history = history + time.asctime() + "." 
##
##    nc_unit = Nio.open_file(nc_file, mode="w",
##                            options=opt, history=history )
##
##    #----------------------------------------------
##    # Create grid dimensions nx and ny, plus time
##    #----------------------------------------------
##    nc_unit.create_dimension("nx", nx)
##    nc_unit.create_dimension("ny", ny)
##    nc_unit.create_dimension("time", None)   # (unlimited dimension)
##    
##    #--------------------------------
##    # Create a variable in the file
##    #----------------------------------
##    # Returns "var" as a PyNIO object
##    #----------------------------------
##    dtype = "d"    # (double, Float64)
##    # dtype = "f"  # (float,  Float32)
##    # dtype = "l"  # (long,   Int64)
##    # dtype = "i"  # (int,    Int32)
##    # dtype = "h"  # (short,  Int16)
##    # dtype = "b"  # (byte,   Int8)
##    # dtype = "S1" # (char)
##
##    var = nc_unit.create_variable(var_name, dtype, ("time", "ny", "nx"))
##    ## var = nc_unit.create_variable(var_name, dtype, ("time", "nx", "ny"))
##
##    #-------------------------------------------
##    # Create a separate, scalar "time stamp" ?
##    #-------------------------------------------
##    # t = nc_unit.create_variable("time", dtype, ("time"))
##    
##    #----------------------------------
##    # Specify a "nodata" fill value ?
##    #----------------------------------
##    var._FillValue = -9999.0    ## Does this jive with Prefill above ??
##    
##    #------------------------------------
##    # Create attributes of the variable
##    #------------------------------------
##    nc_unit.variables[var_name].long_name = long_name
##    nc_unit.variables[var_name].units     = units_name
##
##    return nc_unit
##
###   open_new_netcdf_file()
###-------------------------------------------------------------------
##def add_grid_to_netcdf_file(nc_unit, var_name, grid, time_index):
##
##    #---------------------------------
##    # Assign a value to the variable
##    #-------------------------------------------
##    # This syntax works for scalars and grids
##    #-------------------------------------------
##    # nc_unit.variables[var_name].assign_value( grid )
##
##    var = nc_unit.variables[var_name]
##    var[time_index] = grid
##
##    #-----------------------------------------
##    # Print a summary of the file's contents
##    #-----------------------------------------
##    # print nc_unit
##
###   add_grid_to_netcdf_file()
###-------------------------------------------------------------------
##def get_grid_from_netcdf_file(nc_unit, var_name, time_index):
##
##    var = nc_unit.variables[ var_name ]
##    return var[ time_index ]
##    
###   get_grid_from_netcdf_file()
###-------------------------------------------------------------------
##def close_netcdf_file(nc_unit):
##
##    nc_unit.close()
##
###   close_netcdf_file()
###-------------------------------------------------------------------
##def write_grid_to_netcdf(nc_file, grid, var_name,
##                         long_name, units_name):
##
##    #----------------------------------------------------
##    # Notes:  See: http://www.pyngl.ucar.edu/Nio.shtml.
##    #         Make sure to "import Nio" at the top.
##    #----------------------------------------------------
##    # To retrieve information from Nio object:
##    #    var      = f.variables[var_name]
##    #    type     = var.typecode()
##    #    numDims  = var.rank
##    #    dimSizes = var.shape
##    #    dimNames = var.dimensions
##    #----------------------------------------------------
##    # 8/11/09, but unresolved "gfortran library issue
##    #    on Mac OS X.
##    #----------------------------------------------------
##    # 10/6/09.  On beach, use Python 2.6 at:
##    #   /usr/local/python/bin/python
##    #----------------------------------------------------
##    try:
##        import Nio  # (a module in the PyNIO package)
##    except:
##        python_version = sys.version[:3]
##        print 'SORRY, Cannot write netCDF files because'
##        print 'the "Nio" module cannot be imported.'
##        print 'Note that PyNIO is only installed for'
##        print 'Python version 2.6 on "beach".'
##        print 'The current Python version is:' + python_version
##        print ' '
##        return
##
##    #-------------------------------------
##    # Open a new netCDF file for writing
##    #-------------------------------------
##    opt = Nio.options()
##    opt.PreFill = False
##    f = Nio.open_file(nc_file, mode="w",
##                      options=opt,
##                      history="Created by TopoFlow.")
##
##    #-----------------------------------
##    # Create grid dimensions nx and ny
##    #-----------------------------------
##    ny = numpy.size(grid, 0)
##    nx = numpy.size(grid, 1)
##    f.create_dimension('nx', nx)
##    f.create_dimension('ny', ny)
##
##    #--------------------------------
##    # Create a variable in the file
##    #--------------------------------
##    # Returns "var" as an Nio object
##    #--------------------------------
##    # var = f.create_variable(var_name, 'f', ('nx','ny'))
##    var = f.create_variable(var_name, 'd', ('nx','ny'))
##    
##    #------------------------------------
##    # Create attributes of the variable
##    #------------------------------------
##    f.variables[var_name].long_name = long_name
##    f.variables[var_name].units     = units_name
##
##    #-------------------------------------------
##    # Assign a value to the variable and close
##    #-------------------------------------------
##    # This syntax works for scalars and grids
##    #-------------------------------------------
##    f.variables[var_name].assign_value( grid )
##    f.close()
##    
###   write_grid_to_netcdf()
###-------------------------------------------------------------------  
##def write_pixel_file_header(file_unit, outlet_IDs,
##                            ####  nx,
##                            var_name=None):
##
##    #------------------------------------------------------
##    # Notes:  This is currently hardwired to have columns
##    #         that are each 15 characters wide.
##    #------------------------------------------------------
##    if (var_name is None): var_name = 'F'
##
##    #---------------------------------------------------------
##    # NB!  Allow "outlet_IDs" to be passed as either a tuple
##    #      (numpy.where style indices) or as a 1D array
##    #      (calendar style indices).
##    #---------------------------------------------------------
##    if ('tuple' in str(type(outlet_IDs))):
##        rows = outlet_IDs[0]
##        cols = outlet_IDs[1]
##        print 'In model_output.write_pixel_file_header():'
##        print 'outlet_rows =', rows
##        print 'outlet_cols =', cols
##        print ' '
####    else:
####        rows = (outlet_IDs / nx)
####        cols = (outlet_IDs % nx)
##        
##    n_outlets = size(rows)
##    col_width = 15
##    
##    #--------------------------
##    # Create the heading text
##    #--------------------------
##    file_unit.write( 'Time [min]'.rjust(col_width) )
##    for k in xrange(n_outlets):
##        str1 = str(rows[k])
##        str2 = str(cols[k])
##        str3 = var_name + '[' + str1 + ',' + str2 + ']'
##        file_unit.write( str3.rjust(col_width) ) 
##    file_unit.write("\n")
##    
##    #------------------------------------
##    # Draw a horizontal line of hyphens
##    #------------------------------------
##    width = (n_outlets + 1) * col_width
##    hline = ''.ljust(width, '-')
##    file_unit.write(hline + "\n")
##    
###   write_pixel_file_header()
###-------------------------------------------------------------------
##def write_pixel_file_line(file_unit, time_min, var, outlet_IDs):
##
##    #---------------------------------------------------------
##    # NB!  Allow "outlet_IDs" to be passed as either a tuple
##    #      (numpy.where style indices) or as a 1D array
##    #      (calendar style indices).
##    #---------------------------------------------------------
##
##    #--------------------------------------------------------
##    # Assume that outlet_IDs are tuples as with numpy.where
##    # Note: ('tuple' in str(type(outlet_IDs))) is True.
##    #--------------------------------------------------------
##    rows = outlet_IDs[0]
##    cols = outlet_IDs[1]
##
##    #------------------------------------------
##    # Assume that outlet_IDs is a 1D array of
##    # long-integer "calendar-style" indices
##    #------------------------------------------
####        rows = (outlet_IDs / nx)
####        cols = (outlet_IDs % nx)
##        
##    n_outlets = size(rows)
##    col_width = 15
##
##    #--------------------------------------------------------
##    # If "var" is a grid, subscript with outlet_ID to get a 1D
##    # array of values.  If "var" is scalar, return a vector
##    # with the scalar value repeated once for each "outlet".
##    #--------------------------------------------------------    
##    vals = Pixel_Var(var, outlet_IDs)
##    
##    tstr = ('%15.7f' % time_min)
##    file_unit.write( tstr.rjust(col_width)  )
##    for k in xrange(n_outlets):
##        vstr = ('%15.7f' % vals[k])   ## (same as v. 1.5b)
##        file_unit.write( vstr.rjust(col_width) ) 
##    file_unit.write("\n")
##    
###   write_pixel_file_line()   
###-------------------------------------------------------------------
