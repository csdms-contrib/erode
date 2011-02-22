
# S.D. Peckham
# October 14, 2009

import os
import os.path
import sys

import numpy
import bov_files
import rti_files

#-------------------------------------------------------------------
#
#   unit_test()
#
#   class rts_file():
#
#       open_file()
#       open_new_file()
#       add_grid()
#       get_grid()
#       close_file()
#       close()
#       --------------------
#       byte_swap_needed()
#       number_of_grids()
#
#-------------------------------------------------------------------
def unit_test(nx=4, ny=5, n_grids=6, VERBOSE=False,
              file_name="TEST_FILE.rts"):

    print 'Running unit_test()...'

    #------------------------------------
    # Make instance of rts_file() class
    #------------------------------------
    rts = rts_file()  
    dx = 100
    dy = 100

    #---------------------------------
    # These are unused for RTS files
    #---------------------------------
##    grid_name  = "depth"
##    long_name  = "depth of water"
##    units_name = "meters"

    info = rti_files.make_info( file_name, nx, ny, dx, dy )
    OK = rts.open_new_file( file_name, info )

    if not(OK):
        print 'ERROR during open_new_file().'
        return
    
    grid = numpy.arange(nx * ny, dtype='Float32')
    grid = grid.reshape( (ny, nx) )
    
    #----------------------------------
    # Add some test grids to the file
    #----------------------------------
    for time_index in xrange(n_grids):
        rts.add_grid( grid )                         
        grid = (grid + 1)
        
    rts.close_file()
    print 'Finished writing file: ' + file_name
    print ' '

    #---------------------------------------------
    # Re-open the file and read grids one-by-one 
    #---------------------------------------------
    OK = rts.open_file( file_name )
    if not(OK): return
    n_grids = rts.number_of_grids()
    print 'Reading grids from RTS file: '
    print 'rts.number_of_grids()  =', n_grids
    print 'rts.byte_swap_needed() =', rts.byte_swap_needed()
    print ' '
    for time_index in xrange(n_grids):
        grid = rts.get_grid( time_index )
        print 'grid[' + str(time_index) + '] = '
        print grid
        print '-----------------------------------------------'

    #----------------------------
    # Go back and read 2nd grid
    #----------------------------
    grid = rts.get_grid( 1 )
    print ' '
    print 'Reading second grid again...'
    print 'Second grid ='
    print grid
    print '-----------------------------------------------'
    rts.close_file()
    print 'Finished reading file: ' + file_name
    print ' '

    #---------------------------------------
    # Re-open the file and change one grid
    #---------------------------------------
    print 'Updating RTS file:', file_name
    grid = numpy.ones( (ny, nx), dtype='Float32' )
    OK = rts.open_file( file_name, UPDATE=True )
    if not(OK): return
    rts.add_grid( grid, time_index=0 )
    rts.close_file()
    print 'Finished updating RTS file.'
    print ' '
    
    #---------------------------------------------
    # Re-open the file and read grids one-by-one 
    #---------------------------------------------
    OK = rts.open_file( file_name )
    if not(OK): return
    n_grids = rts.number_of_grids()
    print 'Reading grids from RTS file: '
    print 'rts.number_of_grids()  =', n_grids
    print 'rts.byte_swap_needed() =', rts.byte_swap_needed()
    print ' '
    for time_index in xrange(n_grids):
        grid = rts.get_grid( time_index )
        print 'grid[' + str(time_index) + '] = '
        print grid
        print '-----------------------------------------------'
    rts.close_file()
    print 'Finished reading file: ' + file_name
    print ' '    
    
#   unit_test()
#-------------------------------------------------------------------
class rts_file():

    #----------------------------------------------------------
    def open_file(self, file_name, UPDATE=False):

        info = rti_files.read_info( file_name )
        if (info == -1): return

        #----------------------
        # Store info in state
        #----------------------
        self.info = info
        self.nx   = info.ncols
        self.ny   = info.nrows
        self.dx   = info.xres
        self.dy   = info.yres

        BPE = rti_files.get_bpe( info.data_type )
        self.grid_size   = (self.nx * self.ny * BPE)
        self.SWAP_ENDIAN = self.byte_swap_needed()
        self.file_name   = file_name
        self.time_index  = 0
        
        #-----------------------------------
        # Open file to read only or update
        #-----------------------------------        
        try:
            if (UPDATE):
                rts_unit = open(file_name, 'rb+')
                self.rts_unit = rts_unit
            else:
                rts_unit = open(file_name, 'rb')
                self.rts_unit = rts_unit
            ### return rts_unit
            return True
        except:
            print 'ERROR during rts.open_file().'
            return False
    
    #   open_file()
    #----------------------------------------------------------
    def check_and_store_info(self, file_name, info=None,
                             var_name='UNKNOWN',
                             dtype='float32',
                             MAKE_RTI=True, MAKE_BOV=False):

        #-----------------------------------------------------
        # Note: This object (self) may be new or it may have
        #       been used previously.  In the latter case,
        #       "info" should still be available in "self".
        #       We only need info if MAKE_RTI or MAKE_BOV.
        #-----------------------------------------------------
        self.format     = 'RTS'
        self.file_name  = file_name
        self.time_index = 0  # (need here for RTS files)
        if not(MAKE_RTI or MAKE_BOV): return

        #---------------------------------
        # Was "info" argument provided ?
        #---------------------------------
        NEW_INFO = True
        if (info == None):
            try:
                info = self.info
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
                    ## print 'Reading info from: ' + RTI_file
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
            
##        #---------------------------------
##        # Was "info" argument provided ?
##        #---------------------------------
##        if (info != None):
##            #------------------------------
##            # Save info to a new RTI file
##            #------------------------------
##            prefix   = rti_files.get_file_prefix( file_name )
##            RTI_file = (prefix + '.rti')
##            rti_files.write_info( RTI_file, info )          
##
##        else:
##            #------------------------------------------
##            # Try to find RTI file to copy info from.
##            # Don't create a new RTI file.
##            #------------------------------------------
##            RTI_file = rti_files.try_to_find_rti_file( file_name )
##            if (RTI_file != 'none'):
##                info = rti_files.read_info( RTI_file )
##                info.file_name = file_name
##                info.data_type = rti_files.get_rti_data_type( dtype )
##            else:
##                print 'ERROR during open_new_file():'
##                print '   Could not find RTI file and "info"'
##                print '   argument was not provided.'
##                print ' '
##                return

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
            bov_files.write_info_as_bov( file_name, info, var_name)
                                         ###  time )
        
    #   check_and_store_info()
    #----------------------------------------------------------
    def open_new_file(self, file_name, info=None,
                      var_name='UNKNOWN', dtype='float32',
                      VERBOSE=False,
                      MAKE_RTI=True, MAKE_BOV=False):

        #---------------------------------------
        # Check and store the grid information
        #---------------------------------------
        self.check_and_store_info( file_name, info, var_name,
                                   dtype, MAKE_RTI, MAKE_BOV )
        
        #------------------------------------
        # Try to open new RTS file to write
        #------------------------------------
        try:            
            if (VERBOSE):
                print 'Preparing to write new RTS file:'
                print '   ' + file_name  
            rts_unit = open(file_name, 'wb')
            self.rts_unit = rts_unit
            return True
        except:
            return False
        
    #   open_new_file()
    #----------------------------------------------------------
    def add_grid(self, grid, time_index=-1):

        #---------------------------------------------
        # Can use time_index to move file pointer
        # and overwrite an existing grid vs. append.
        #---------------------------------------------
        if (time_index >= 0):
            offset = (time_index * self.grid_size)
            self.rts_unit.seek( offset )

        #--------------------------------
        # Swap byte order, if necessary
        #--------------------------------
        if (self.info.SWAP_ENDIAN): grid.byteswap(True)

        #-------------------------------
        # Write grid as binary to file
        #-------------------------------
        grid.tofile( self.rts_unit )
        self.time_index += 1
        
    #   add_grid()
    #----------------------------------------------------------
    def get_grid(self, time_index, dtype='float32'):

        #-----------------------------------------------
        # Compute offset from time_index and grid_size
        #-----------------------------------------------
        n_values = self.nx * self.ny
        offset   = (time_index * self.grid_size)
        self.rts_unit.seek( offset )

        grid = numpy.fromfile( self.rts_unit, count=n_values,
                               dtype=dtype )
        grid = grid.reshape( self.ny, self.nx )

        #--------------------------------
        # Swap byte order, if necessary
        #--------------------------------
        if (self.info.SWAP_ENDIAN): grid.byteswap(True)
            
        return grid
    
    #   get_grid()
    #-------------------------------------------------------------------
    def close_file(self):

        self.rts_unit.close()

    #   close_file()
    #-------------------------------------------------------------------
    def close(self):

        self.rts_unit.close()

    #   close()
    #-------------------------------------------------------------------
    def byte_swap_needed(self):

        machine_byte_order = rti_files.get_rti_byte_order()
        SWAP =  (machine_byte_order != self.info.byte_order)
        
        return SWAP
    
    #   byte_swap_needed()
    #-------------------------------------------------------------------
    def number_of_grids(self):

        file_size = os.path.getsize( self.file_name )
        n_grids   = (file_size / self.grid_size)

        # self.file_size = file_size
        # self.n_grids   = n_grids
        
        return n_grids
        
    #   number_of_grids()
    #-------------------------------------------------------------------

    
