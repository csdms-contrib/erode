
# S.D. Peckham
# October 8, 2009

from numpy import *

import os
import sys
import time
import Nio    # (a module in the PyNIO package) 

#-------------------------------------------------------------------

#   Unit_Test1()
#   open_new_netcdf_file()
#   add_grid_to_netcdf_file()
#   get_grid_from_netcdf_file()   # (assumed already open)
#   close_netcdf_file()

#-------------------------------------------------------------------
def Unit_Test1():

    import Nio    # (a module in the PyNIO package) 

    print 'Running Unit_Test1()...'
    nc_file = "TEST_FILE.nc"
    nx = 4
    ny = 4
    var_name   = "depth"
    long_name  = "depth of water"
    units_name = "meters"
    n_grids    = 5
 
    nc_unit = open_new_netcdf_file( nc_file, nx, ny, var_name, 
                                    long_name, units_name )

    grid = arange(nx * ny, dtype='Float32').reshape((ny,nx))

    for time_index in xrange(n_grids):
        add_grid_to_netcdf_file(nc_unit, var_name, grid, time_index)                         
        grid = (grid + 1)
        print nc_unit  # (print a summary)

    close_netcdf_file( nc_unit )
    ## nc_unit.close()
    print 'Finished writing file: ' + nc_file
    print ' '

    #------------------------------------------
    # Open the file and read grids one-by-one 
    #------------------------------------------
    nc_unit = Nio.open_file(nc_file, mode="r")
    for time_index in xrange(n_grids):
        grid = get_grid_from_netcdf_file(nc_unit, var_name, time_index)
        print 'grid[' + str(time_index) + '] = '
        print grid
        print '-----------------------------------------------'
    close_netcdf_file( nc_unit )

#   Unit_Test1()
#-------------------------------------------------------------------
def open_new_netcdf_file(nc_file, nx, ny,
                         var_name, long_name, units_name):

    try:
        import Nio
        print 'Imported Nio version: ' + Nio.__version__
    except:
        python_version = sys.version[:3]
        print 'SORRY, Cannot write netCDF files because'
        print 'the "Nio" package cannot be imported.'
        print 'Note that "PyNIO" is only installed for'
        print 'Python version 2.6 on "beach".'
        print 'The current Python version is:' + python_version
        print ' '
        return

    #----------------------------
    # Does file already exist ?
    #--------------------------------------------
    # Overwrite not allowed; must remove first.
    #--------------------------------------------
    if (os.path.exists( nc_file )):
        print 'Deleting existing file: ' + nc_file
        os.remove( nc_file )

    #----------------------------
    # Does file already exist ?
    #----------------------------
##    if (os.path.exists(nc_file)):
##        print 'WARNING: Found a file named:'
##        print '  ' + nc_file
##        print 'in the current working directory.'
##        print 'Are you sure you want to overwrite it?'
##        ## os.remove(nc_file)   # (delete existing file)
##        return
        
    #-------------------------------------
    # Open a new netCDF file for writing
    #-------------------------------------
    # Sample output from time.asctime():
    #     "Thu Oct  8 17:10:18 2009"
    #-------------------------------------
    opt = Nio.options()
    opt.PreFill = False            # (for efficiency)
    opt.HeaderReserveSpace = 4000  # (4000 bytes, for efficiency)
    history = "Created by TopoFlow 3.0 on "
    history = history + time.asctime() + "." 

    nc_unit = Nio.open_file(nc_file, mode="w",
                            options=opt, history=history )

    #----------------------------------------------
    # Create grid dimensions nx and ny, plus time
    #----------------------------------------------
    nc_unit.create_dimension("nx", nx)
    nc_unit.create_dimension("ny", ny)
    nc_unit.create_dimension("time", None)   # (unlimited dimension)
    
    #--------------------------------
    # Create a variable in the file
    #----------------------------------
    # Returns "var" as a PyNIO object
    #----------------------------------
    dtype = "d"    # (double, Float64)
    # dtype = "f"  # (float,  Float32)
    # dtype = "l"  # (long,   Int64)
    # dtype = "i"  # (int,    Int32)
    # dtype = "h"  # (short,  Int16)
    # dtype = "b"  # (byte,   Int8)
    # dtype = "S1" # (char)

    var = nc_unit.create_variable(var_name, dtype, ("time", "ny", "nx"))
    ## var = nc_unit.create_variable(var_name, dtype, ("time", "nx", "ny"))

    #-------------------------------------------
    # Create a separate, scalar "time stamp" ?
    #-------------------------------------------
    # t = nc_unit.create_variable("time", dtype, ("time"))
    
    #----------------------------------
    # Specify a "nodata" fill value ?
    #----------------------------------
    var._FillValue = -9999.0    ## Does this jive with Prefill above ??
    
    #------------------------------------
    # Create attributes of the variable
    #------------------------------------
    nc_unit.variables[var_name].long_name = long_name
    nc_unit.variables[var_name].units     = units_name

    return nc_unit

#   open_new_netcdf_file()
#-------------------------------------------------------------------
def add_grid_to_netcdf_file(nc_unit, var_name, grid, time_index):

    #---------------------------------
    # Assign a value to the variable
    #-------------------------------------------
    # This syntax works for scalars and grids
    #-------------------------------------------
    # nc_unit.variables[var_name].assign_value( grid )

    var = nc_unit.variables[var_name]
    var[time_index] = grid

    #-----------------------------------------
    # Print a summary of the file's contents
    #-----------------------------------------
    # print nc_unit

#   add_grid_to_netcdf_file()
#-------------------------------------------------------------------
def get_grid_from_netcdf_file(nc_unit, var_name, time_index):

    var = nc_unit.variables[ var_name ]
    return var[ time_index ]
    
#   get_grid_from_netcdf_file()
#-------------------------------------------------------------------
def close_netcdf_file(nc_unit):

    nc_unit.close()

#   close_netcdf_file()
#-------------------------------------------------------------------
