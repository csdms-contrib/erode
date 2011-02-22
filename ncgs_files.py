
# S.D. Peckham
# October 13, 2009
# December 2, 2009 (updated open_new_file to use "info")

import os
import sys
import time

import numpy
import bov_files
import rti_files

# import Nio    # (a module in the PyNIO package) 

#-------------------------------------------------------------------

#   unit_test()
#   save_ncgs_frame()  #### (12/7/09) ####

#   class ncgs_file():

#       import_nio()
#       open_file()
#       check_and_store_info()   # (12/2/09)
#       open_new_file()
#       add_grid()
#       get_grid()
#       close_file()
#       close()

#-------------------------------------------------------------------
def unit_test(nx=4, ny=5, n_grids=6, VERBOSE=False,
              file_name="TEST_FILE.nc"):

    print 'Running unit_test()...'

    #-------------------------------------
    # Make instance of ncgs_file() class
    #-------------------------------------
    ncgs = ncgs_file()
    dx = 100
    dy = 100
    grid_name = "depth"

    info = rti_files.make_info( file_name, nx, ny, dx, dy )
    OK = ncgs.open_new_file( file_name, info,
                             dtype='float32',
                             grid_name=grid_name,
                             long_name="depth of water",
                             units_name="meters",
                             comment="Created by TopoFlow 3.0.")
    if not(OK):
        print 'ERROR during open_new_file().'
        return
    
    grid = numpy.arange(nx * ny, dtype='Float32')
    grid = grid.reshape( (ny, nx) )

    #----------------------------------
    # Add some test grids to the file
    #----------------------------------
    print 'Writing grids to NCSG file...'
    for time_index in xrange(n_grids):
        ncgs.add_grid( grid, grid_name )  
        ## ncgs.add_grid( grid, grid_name, time_index )                         
        grid = (grid + 1)
        
    if (VERBOSE):
        print self.ncgs_unit  # (print a summary)

    ncgs.close_file()
    print 'Finished writing NCGS file: ' + file_name
    print ' '

    #---------------------------------------------
    # Re-open the file and read grids one-by-one 
    #---------------------------------------------
    OK = ncgs.open_file( file_name )
    if not(OK): return
    print 'Reading grids from NCGS file: '
    
    for time_index in xrange(n_grids):
        grid = ncgs.get_grid(grid_name, time_index)
        print 'grid[' + str(time_index) + '] = '
        print grid
        print '-----------------------------------------------'
    ncgs.close_file()    
    print 'Finished reading NCGS file: ' + file_name
    print ' '
    
#   unit_test()
#-------------------------------------------------------------------
def save_ncgs_frame(ncgs_file_name=None, rtg_file_name=None):

    ncgs = ncgs_file()
    OK = ncgs.open_file( ncgs_file_name )
    if not(OK): return

    grid_name  = 'H'
    time_index = 200
    grid = ncgs.get_grid( grid_name, time_index )
    ncgs.close()
    
    grid = numpy.array( grid )
    print 'min(grid), max(grid) =', grid.min(), grid.max()

    rtg_unit = open( rtg_file_name, 'wb' )
    grid.tofile( unit )
    rtg_unit.close()

#   save_ncgs_frame()
#-------------------------------------------------------------------
class ncgs_file():

    #----------------------------------------------------------
    # Note:  ncgs = NetCDF Grid Stack (used by CSDMS)
    #----------------------------------------------------------
    def import_nio(self):

        try:
            import Nio  # (a module in the PyNIO package) 
            print 'Imported Nio version: ' + Nio.__version__
            return Nio
        except:
            python_version = sys.version[:3]
            print ' '
            print 'SORRY, Cannot write netCDF files because'
            print 'the "Nio" package cannot be imported.'
            print ' '
            if (python_version != '2.6'):
                print 'Note that "PyNIO" is only installed for'
                print 'Python version 2.6 on "beach".'
                print 'The current Python version is:', python_version
                print ' '
            return False
        
    #   import_nio()
    #----------------------------------------------------------
    def open_file(self, file_name):

        #--------------------------------------------------
        # Try to import the Nio module from PyNIO package
        #--------------------------------------------------
        Nio = self.import_nio()
        if not(Nio): return
        
        #-------------------------
        # Open file to read only
        #-------------------------
        try:
            ncgs_unit = Nio.open_file(file_name, mode="r")
            self.ncgs_unit = ncgs_unit
            ### return ncgs_unit
            return True
        except:
            return False
    
    #   open_file()
    #----------------------------------------------------------
    def check_and_store_info(self, file_name, info=None,
                             grid_name='UNKNOWN',
                             dtype='float32',
                             MAKE_RTI=True, MAKE_BOV=False):

        #-----------------------------------------------------
        # Note: This object (self) may be new or it may have
        #       been used previously.  In the latter case,
        #       "info" should still be available in "self".
        #       We only need info if MAKE_RTI or MAKE_BOV.
        #-----------------------------------------------------
        self.format     = 'NCGS'
        self.file_name  = file_name
        self.time_index = 0
        self.grid_name  = grid_name
        
        #-----------------------------------------------------
        # This was used by rts_files.check_and_store_info()
        # but is not appropriate here because we need to
        # know nx, ny, dx and dy for the netCDF file.
        #-----------------------------------------------------        
        ### if not(MAKE_RTI or MAKE_BOV): return

        #---------------------------------
        # Was "info" argument provided ?
        #---------------------------------
        NEW_INFO = True
        if (info == None):
            try:
                info    = self.info
                self.nx = info.ncols  ###
                self.ny = info.nrows
                NEW_INFO = False
                ## print 'Found info in state.'
            except:
                #------------------------------------------
                # Try to find RTI file to copy info from.
                # Don't create a new RTI file.
                #------------------------------------------
                RTI_file = rti_files.try_to_find_rti_file( file_name )
                if (RTI_file != 'none'):
                    info = rti_files.read_info( RTI_file )
                    print 'Reading info from: ' + RTI_file
                else:
                    print 'ERROR during open_new_file():'
                    print '   Could not find RTI file and "info"'
                    print '   argument was not provided.'
                    print ' '
                    return

        #-----------------------------
        # Update "info" as necessary
        #-----------------------------
        info.grid_file   = file_name
        info.data_type   = rti_files.get_rti_data_type( dtype )
        info.data_source = 'TopoFlow 3.0'
        info.gmin        = -9999.0
        info.gmax        = -9999.0
        
        #---------------------------------------
        # If new "info" was provided, store it
        #---------------------------------------
        if (NEW_INFO):
            self.info = info
            self.nx   = info.ncols
            self.ny   = info.nrows           
            ## print 'Stored new info in state.'

        #-------------------
        # Write RTI file ?
        #-------------------
        if (MAKE_RTI):
            prefix   = rti_files.get_file_prefix( file_name )
            RTI_file = (prefix + '.rti')
            rti_files.write_info( RTI_file, info )
            # print 'Wrote grid info to: ' + RTI_file   ######
            
        #-------------------
        # Write BOV file ?
        #-------------------
        if (MAKE_BOV):
            bov_files.write_info_as_bov( file_name, info, grid_name)
                                         ###  time )
        
    #   check_and_store_info()
    #----------------------------------------------------------
    def open_new_file(self, file_name, info=None,
                      grid_name='X',
                      long_name=None,
                      units_name='None',
                      dtype='float32',
                      comment='',
                      MAKE_RTI=True, MAKE_BOV=False):

        #----------------------------------------
        # Possible settings for "nio_type_code"
        #-------------------------------------------
        # nio_type_code = "d"  # (double, Float64)
        # nio_type_code = "f"  # (float,  Float32)
        # nio_type_code = "l"  # (long,   Int64)
        # nio_type_code = "i"  # (int,    Int32)
        # nio_type_code = "h"  # (short,  Int16)
        # nio_type_code = "b"  # (byte,   Int8)
        # nio_type_code = "S1" # (char)
        #-------------------------------------------
        nio_type_map = {'float64':'d', 'float32':'f',
                        'int64':'l', 'int32':'i',
                        'int16':'s', 'int8':'b',
                        'S|100':'S1'}  # (check last entry)
        nio_type_code = nio_type_map[ dtype.lower() ]
            
        #--------------------------------------------------
        # Try to import the Nio module from PyNIO package
        #--------------------------------------------------
        Nio = self.import_nio()
        if not(Nio): return False

        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.check_and_store_info( file_name, info, grid_name,
                                   dtype, MAKE_RTI, MAKE_BOV )
        if (long_name == None): long_name = grid_name
        self.long_name  = long_name
        self.units_name = units_name
        self.dtype      = dtype
        self.nio_type_code = nio_type_code   ######
        
        #----------------------------
        # Does file already exist ?
        #--------------------------------------------
        # Overwrite not allowed; must remove first.
        #--------------------------------------------
        if (os.path.exists( file_name )):
            print 'Deleting existing file: ' + file_name
            os.remove( file_name )

        #----------------------------
        # Does file already exist ?
        #----------------------------
    ##    if (os.path.exists( file_name )):
    ##        print 'WARNING: Found a file named:'
    ##        print '  ' + file_name
    ##        print 'in the current working directory.'
    ##        print 'Are you sure you want to overwrite it?'
    ##        ## os.remove( file_name )   # (delete existing file)
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
        history = "Created using PyNIO " + Nio.__version__ + " on "
        history = history + time.asctime() + ". " 
        history = history + comment

        try:
            ncgs_unit = Nio.open_file(file_name, mode="w",
                                      options=opt, history=history )
            OK = True
        except:
            OK = False
            return OK

##        print 'nx =', self.info.ncols
##        print 'ny =', self.info.nrows
##        print 'dx =', self.info.xres
##        print 'dy =', self.info.yres
##        print ' '
        
        #----------------------------------------------
        # Create grid dimensions nx and ny, plus time
        #----------------------------------------------
        # Without usint "int()" here, we get this:
        #     TypeError: size must be None or integer
        #----------------------------------------------
        ncgs_unit.create_dimension("nx", int(self.info.ncols))
        ncgs_unit.create_dimension("ny", int(self.info.nrows))
        ncgs_unit.create_dimension("time", None)   # (unlimited dimension)
        
        #--------------------------------
        # Create a variable in the file
        #----------------------------------
        # Returns "var" as a PyNIO object
        #----------------------------------
        var = ncgs_unit.create_variable(grid_name, nio_type_code,
                                        ("time", "ny", "nx"))
        ## var = nc_unit.create_variable(grid_name, nio_type_code,
        ##            ("time", "nx", "ny"))

        #-------------------------------------------
        # Create a separate, scalar "time stamp" ?
        #-------------------------------------------
        # t = nc_unit.create_variable("time", nio_type_code, ("time"))
        
        #----------------------------------
        # Specify a "nodata" fill value ?
        #----------------------------------
        var._FillValue = -9999.0    ## Does this jive with Prefill above ??
        
        #------------------------------------
        # Create attributes of the variable
        #------------------------------------
        ncgs_unit.variables[grid_name].long_name = long_name
        ncgs_unit.variables[grid_name].units     = units_name
        ncgs_unit.variables[grid_name].dx        = self.info.xres
        ncgs_unit.variables[grid_name].dy        = self.info.yres  ### (12/2/09)
##        ncgs_unit.variables[grid_name].dx        = dx
##        ncgs_unit.variables[grid_name].dy        = dy  ### (10/15/09)
        ncgs_unit.variables[grid_name].y_south_edge = self.info.y_south_edge
        ncgs_unit.variables[grid_name].y_north_edge = self.info.y_north_edge
        ncgs_unit.variables[grid_name].x_west_edge  = self.info.x_west_edge
        ncgs_unit.variables[grid_name].x_east_edge  = self.info.x_east_edge        
        
        self.ncgs_unit = ncgs_unit
        #### return ncgs_unit
        return OK
    
    #   open_new_file()
    #----------------------------------------------------------
    def add_grid(self, grid, grid_name, time_index=-1):

        #---------------------------------
        # Assign a value to the variable
        #-------------------------------------------
        # This syntax works for scalars and grids
        #-------------------------------------------
        # nc_unit.variables[var_name].assign_value( grid )


        #-------------------------------------
        # Can use time_index to overwrite an
        # existing grid vs. simple append.
        #-------------------------------------
        if (time_index == -1):
            time_index = self.time_index

        #---------------------------------------
        # Write a grid to existing netCDF file
        #---------------------------------------
        var = self.ncgs_unit.variables[ grid_name ]
        var[ time_index ] = grid
        self.time_index += 1
        
        #-------------------------------------------------
        # 12/2/09:  netCDF is supposed to take care of
        # byteorder transparently.  However, we need to
        # make sure we don't byteswap in the function
        # "model_output.save_as_grid_to_file()" when the
        # output format is netCDF.
        #-------------------------------------------------        
##        if (sys.byteorder == 'big'):
##            var[time_index] = grid
##        else:
##            grid2 = grid.copy()
##            var[time_index] = grid2.byteswap() 
##        self.time_index += 1
        
    #   add_grid()
    #----------------------------------------------------------
    def get_grid(self, grid_name, time_index):

        var = self.ncgs_unit.variables[ grid_name ]
        return var[ time_index ]
        
    #   get_grid()
    #-------------------------------------------------------------------
    def close_file(self):

        self.ncgs_unit.close()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        self.ncgs_unit.close()

    #   close()
    #-------------------------------------------------------------------
    
