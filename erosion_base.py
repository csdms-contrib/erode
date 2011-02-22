
## Copyright (c) 2001-2009, Scott D. Peckham
## January 2009  (converted from IDL)
## September, October, November 2009

#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class erosion_component
#      __init__()
#-----------------------------
#      initialize()
#      update()
#      finalize()
#      read_config_file() 
#-----------------------------
#      get_cca_ports()           # (in CSDMS_base.py)
#      embed_child_components()
#      add_child_ports()
#      initialize_ports()
#      release_cca_ports()       # (in CSDMS_base.py)
#----------------------------------
#      initialize_d8_vars()
#      initialize_DEM()
#      create_initial_DEM()
#      initialize_boundary_conditions()
#      initialize_computed_vars()
#----------------------------------
#      update_R()
#      update_R_integral()
#      update_U()
#      update_U_integral()
#      update_base_level()
#      update_DEM_edge_values()    ###
#--------------------------------
#      update_d8_vars()
#      update_slope_grid()
#      update_Q_grid()
#      update_Qs_grid()
#      update_DEM()
#      fill_pits_in_DEM()          ###
#      update_DEM_min_and_max()
#      check_stability()
#      check_finished()
#      check_steady_state()        ### (not written yet)
#--------------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#--------------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#---------------------------
#  THESE ARE NOT USED YET
#---------------------------
#  Get_Timestep()
#  Stable_Timestep()
#
#-----------------------------------------------------------------------

from numpy import *
import numpy

import os.path

import cfg_files as cfg
import CSDMS_base
import d8_base
import fill_pits  #####
import midpoints  #####
import model_input
import rtg_files
import tf_utils

from model_output import *

#-----------------------------------------------------------------------
def unit_test(n_steps=10):

    #-------------------------------------------------
    # NOTE: The Treynor_Iowa DEM has no depressions!
    #-------------------------------------------------  
##    directory   = tf_utils.TF_Test_Directory()
##    data_prefix = tf_utils.TF_Test_Data_Prefix()
##    case_prefix = tf_utils.TF_Test_Case_Prefix()

    directory   = '/Applications/Erode/Data/Test1/'
    ## directory   = '/data/progs/erode/3.0/data/Test1/'
    data_prefix = 'Test1'
    case_prefix = 'Test1'
    
    ec = erosion_component()
    ec.CCA   = False
##    ec.SILENT = False
##    ec.REPORT = True
    
    ec.run_model(directory=directory,
                 data_prefix=data_prefix,
                 case_prefix=case_prefix,
                 n_steps=n_steps)
    
    #-------------------------------------------
    # Call initialize() and call update() once
    #-------------------------------------------
##    print 'STATUS =', c.get_status()
##    c.initialize(directory=directory,
##                 data_prefix=data_prefix,
##                 case_prefix=case_prefix, mode="main")
##    print 'STATUS =', c.get_status()
##    time = float64(0)
##    c.update(time)
##    print 'STATUS =', c.get_status()

#   unit_test()
#-----------------------------------------------------------------------
class erosion_component(CSDMS_base.CSDMS_component):

    def __init__(self):

        # print "Instantiating erosion_base component..."
        self.CCA    = tf_utils.TF_Use_CCA()
        self.DEBUG  = False
        self.SILENT = True
        self.REPORT = False

        #------------------------
        # Define some constants
        #------------------------
        self.secs_per_year = float64(31536000)
        self.mm_per_m      = float64(1000)
        self.dz_tolerance  = float64(1e-6)
        
##        self.dz_tolerance  = float64(1e-5)
##        self.dz_tolerance  = float64(1e-4)
##        self.dz_tolerance  = float64(1e-3)

        self.status = 'created'    # (OpenMI 2.0 convention)

    #   __init__()
    #-------------------------------------------------------------------
    def get_input_items(self):

        items = ['z0', 'R', 'U']
        return items
        ## return numpy.array(items)    # (string array vs. list)

    #   get_input_items()
    #-------------------------------------------------------------------
    def get_output_items(self):

        items = ['dt', 'nx', 'ny', 'dx', 'dy', 'dd', 'da', 'ds', 'dw',
                 'R', 'U', 'codes', 'S', 'A', 'Q', 'Qs', 'dz', 'z',
                 'm', 'n', 'K', 'p', 'BLR', 'dz_max',
                 'time_sec', 'time_min']

        return items
        ## return numpy.array(items)    # (string array vs. list)
    
    #   get_output_items()
    #-------------------------------------------------------------------
    def initialize(self, directory=None,
                   data_prefix=None,
                   case_prefix=None, mode="module",
                   SILENT=True):

        #---------------------------------------------------
        # When a user clicks on a component's "run" button,
        # that component's "run_model()" method is called
        # which sets "mode" to "main" and passes the mode
        # to initialize().  If this initialize() method is
        # called by another component (a driver), then the
        # mode has the default setting of "module".  In
        # this case, the directory, etc. is set by the
        # caller, so after "load_user_input()" below, we
        # call "set_directory()" again to override and new
        # directory setting from the user.
        #---------------------------------------------------
        self.mode   = mode
        self.status = 'initializing'  # (OpenMI 2.0 convention)
        if not(SILENT):
            print 'Erosion component: Initializing...'
        self.set_directory(directory, data_prefix, case_prefix)

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.read_config_file() # (Sets nx, ny, dx, dy, etc.)
        self.read_grid_info()   # (stores rti in self, adds "da")
        ### self.FILL_PITS = False
        ############################
        self.USER_SET_VALUES = True     #########################
        if (self.USER_SET_VALUES):        
            #-----------------------------------------
            # This may override values from CFG file
            #-----------------------------------------
            self.load_user_input()
            if (self.mode == "module"):
                self.set_directory(directory, data_prefix, case_prefix)
                
        #---------------------------------------------
        # Create an RTI file with grid info which
        # can be used by this and other components ?
        #---------------------------------------------
        if (self.make_z0_method != 'READ_FILE'):
            DEM_file     = (self.case_prefix + '_2D-z0.rtg')
            RTI_file     = (self.case_prefix + '.rti')
            self.z0_file = DEM_file
            self.rti = rti_files.make_info( DEM_file, self.nx,
                                            self.ny, self.dx, self.dy )
            rti_files.write_info(RTI_file, self.rti)

        #----------------------------------------
        # Has this component been "turned off?"
        #----------------------------------------
        if (self.method == 0):
            self.SAVE_Z_GRIDS  = False    # (It is True by default.)
            self.SAVE_Z_PIXELS = False    # (It is True by default.)
            self.DONE = True
            self.status = 'initialized'  # (OpenMI 2.0 convention) 
            return
                                                  
        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        # Can't move read_input_files() to start of
        # update(), since initial values needed here.
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        #-----------------------
        # Initialize variables
        #-----------------------
        self.initialize_d8_vars()  # (depend on D8 flow grid)
        self.initialize_DEM()
        self.initialize_boundary_conditions()
        self.initialize_computed_vars()
        self.initialize_time_vars( units='years' )
        
        #--------------------------------------------------
        # Make sure self.Q_ts_file is not NULL (12/22/05)
        
        # This is only output file that is set by default
        # and is still NULL if user hasn't opened the
        # output var dialog for the channel process.
        #--------------------------------------------------
        if (self.SAVE_Z_PIXELS and (self.z_ts_file == '')):    
            self.z_ts_file = (case_prefix + '_0D-z.txt')       

        #-----------------------------------------------------
        #  If this component is running in stand-alone mode,
        #  then it must initialize the process modules that
        #  it is going to use.
        #-----------------------------------------------------
        self.initialize_required_components( mode )  #### NOT READY YET ######
        
        #-------------------------------
        # Save these for writing output
        #-------------------------------
        ## self.store_outlet_IDs()   #### NOT READY YET ####
        outlet_rows = numpy.array([0, self.ny/2])
        outlet_cols = numpy.array([0, self.nx/2])
        self.outlet_IDs = (outlet_rows, outlet_cols)
        
        self.open_output_files()
        
        self.status = 'initialized'  # (OpenMI 2.0 convention) 
        
    #   initialize()
    #-------------------------------------------------------------------
    def update(self, time=None, SILENT=None, REPORT=None):

        if (SILENT == None): SILENT=self.SILENT
        if (REPORT == None): REPORT=self.REPORT

        #### if not(SILENT) and (self.time_index == 0):
        if (self.time_index == 0):
            print 'Erosion component: Processing...'
        self.status = 'updating'  # (OpenMI 2.0 convention)

        #---------------------------------
        # Print dz_max to track progress
        #---------------------------------
        if (self.mode == 'main'):
            self.print_time_and_value(self.dz_max, 'dz_max', '[m]',
                                      interval=5.0, PRINT_INDEX=True)
            
        #--------------------------------------
        # Update values from other components
        #--------------------------------------
        self.update_R()
        self.update_R_integral()
        self.update_U()
        self.update_U_integral()
        
        #-------------------------
        # Update computed values
        #-------------------------
        self.update_base_level()
        self.update_DEM_edge_values()   ######
        if (self.FILL_PITS):
            self.fill_pits_in_DEM(SILENT=SILENT)    ############

        #--------------------------------------------
        # Update the D8 flow grid and all vars that
        # depend on it, including D8 area grid.
        #--------------------------------------------
        self.update_d8_vars(SILENT=SILENT, REPORT=REPORT)  #########
        self.update_slope_grid(SILENT=SILENT, REPORT=REPORT)
        self.update_Q_grid(SILENT=SILENT, REPORT=REPORT)
        self.update_Qs_grid(SILENT=SILENT, REPORT=REPORT)
        self.update_DEM(SILENT=SILENT, REPORT=REPORT)
        self.update_DEM_min_and_max()

        ########################################
        # CAN THE UPDATED DEM HAVE PITS ??
        # IF NOT, DON'T CALL FILL_PITS.
        ########################################
        
        #------------------------
        # Check computed values
        #------------------------
        OK = self.check_stability()

        #-------------------------------------------
        # Read from files as needed to update vars 
        #-----------------------------------------------------
        # NB! This is currently not needed for the "erosion
        # process" because values don't change over time and
        # read_input_files() is called by initialize().
        #-----------------------------------------------------
        # if (self.time_index > 0):
        #     self.read_input_files()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        self.write_output_files( time )
        if (OK):
            self.status = 'updated'  # (OpenMI 2.0 convention)
        else:
            self.status = 'failed'
            self.DONE   = True

        #------------------------
        # Update internal clock
        #------------------------
        self.update_time()
        ## print 'time_index =', self.time_index
        
        #-------------------------------------------
        # Check for steady-state condition instead
        #-------------------------------------------
        self.check_finished()   ######################
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI)   
        self.close_input_files()   ##  TopoFlow input "data streams"
        self.close_output_files()

        #-----------------------------
        # Save final DEM to DEM_file        # (write a separate function ###########)
        #-----------------------------
##        DEM_unit = open(self.final_DEM_file, 'wb')
##        if (self.rti.SWAP_ENDIAN):
##            final_DEM = self.DEM.byteswap(True)
##            final_DEM.tofile(DEM_unit)
##        else:
##            self.DEM.tofile( DEM_unit )
##        DEM_unit.close()
    
        self.status = 'finalized'  # (OpenMI)

        #----------------------
        # Print final message
        #----------------------
        print 'Erosion component: Finished.'
        print 'Number of timesteps completed =', self.time_index
        print 'min(z), max(z)   =', self.DEM.min(), self.DEM.max()
        print ' '
        print 'dz_max_vec[0]    =', self.dz_max_vec[0]
        print 'dz_max_vec[nt-1] =', self.dz_max_vec[-1]
        print 'dz_max_vec.min() =', self.dz_max_vec.min()
        print 'dz_max_vec.max() =', self.dz_max_vec.max()
        print ' '
        self.print_run_time('Erosion component')
    
        #---------------------------
        # Release all of the ports
        #----------------------------------------
        # Make this call in "finalize()" method
        # of the component's CCA Imple file
        #----------------------------------------
        # self.release_cca_ports( port_names, d_services )
        
    #   finalize()
    #-------------------------------------------------------------------
    def read_config_file(self):
   
        #------------------------------------------
        # Read parameters from a CFG file that is
        # in the current working directory.
        #------------------------------------------
        print 'Erosion component: Reading config file...'
        self.config_file = (self.directory +
                            self.case_prefix + '_erosion.cfg')
        file_unit = open(self.config_file, 'r')
        
        #-------------------------
        # Skip over header lines
        #-------------------------
        cfg.skip_header( file_unit, n_lines=4 )
            
        #------------------------
        # Read the erosion vars
        #------------------------
        self.method      = cfg.read_value( file_unit, dtype='int16' )
        self.method_name = cfg.read_value( file_unit, dtype='string' )
        #-----------------------------------------------------------------        
        dt_info          = cfg.read_input_option( file_unit )
        # self.dt_type   = dt_info[0]
        self.dt          = dt_info[1]
        #-----------------------------------------------------------
        # Read method to get initial elevations, z0
        # Read m, n, K, p, R, U, BLR, etc.
        #-----------------------------------------------------------
        self.nx      = cfg.read_value( file_unit, dtype='int32' )
        self.ny      = cfg.read_value( file_unit, dtype='int32' )
        self.n_steps = cfg.read_value( file_unit, dtype='int32' )
        self.dx      = cfg.read_value( file_unit, dtype='float32' )
        self.dy      = cfg.read_value( file_unit, dtype='float32' )
        self.m       = cfg.read_value( file_unit, dtype='float32' )
        self.n       = cfg.read_value( file_unit, dtype='float32' )
        self.K       = cfg.read_value( file_unit, dtype='float32' )
        self.p       = cfg.read_value( file_unit, dtype='float32' )
        self.R       = cfg.read_value( file_unit, dtype='float32' )
        self.U       = cfg.read_value( file_unit, dtype='float32' ) # (allow grid or grid_stack later)
        self.BLR     = cfg.read_value( file_unit, dtype='float32' ) # (Base-level Lowering Rate)
        #----------------------------------------------------------------------------
##        print 'method      =', self.method
##        print 'method_name =', self.method_name
##        print 'dt          =', self.dt
##        print 'nx          =', self.nx
##        print 'ny          =', self.ny
##        print 'dx          =', self.dx
##        print 'dy          =', self.dy
##        print 'm           =', self.m
##        print 'n           =', self.n
##        print 'K           =', self.K
##        print 'p           =', self.p
##        print 'R           =', self.R
##        print 'U           =', self.U
##        print 'BLR         =', self.BLR
        #-----------------------------------------------------------------------
        self.make_z0_method = cfg.read_value( file_unit, dtype='string' )
        self.FLAT           = (self.make_z0_method.upper() == 'FLAT')
        self.PLANE          = (self.make_z0_method.upper() == 'PLANE')
        self.CORNER_PLANE   = (self.make_z0_method.upper() == 'CORNER_PLANE')
        self.READ_FILE      = (self.make_z0_method.upper() == 'READ_FILE')
        self.z0_plane_dz_dx = cfg.read_value( file_unit, dtype='float32' )
        self.z0_plane_dz_dy = cfg.read_value( file_unit, dtype='float32' )
        self.z0_plane_S     = numpy.sqrt(self.z0_plane_dz_dx**2 + \
                                         self.z0_plane_dz_dy**2)
        self.z0_file        = cfg.read_value( file_unit, dtype='string' )
        #-----------------------------------------------------------------------
##        print 'make_z0_method =', self.make_z0_method
##        print 'z0_plane_dz_dx =', self.z0_plane_dz_dx
##        print 'z0_plane_dz_dy =', self.z0_plane_dz_dy
##        print 'z0_file        =', self.z0_file
        #----------------------------------------------------------------------------
        self.noise_method   = cfg.read_value( file_unit, dtype='string' )
        self.GAUSSIAN       = (self.noise_method.upper() == 'GAUSSIAN')
        self.MIDPOINTS      = (self.noise_method.upper() == 'MIDPOINTS')
        self.NO_NOISE       = (self.noise_method.upper() == 'NO_NOISE')
        self.noise_scale    = cfg.read_value( file_unit, dtype='float32' )
        self.seed           = cfg.read_value( file_unit, dtype='int32' ) # (e.g. 36421)
        #----------------------------------------------------------------------------
##        print 'noise_method   =', self.noise_method
##        print 'seed           =', self.seed
        #----------------------------------------------------------------------------
        self.BC_method      = cfg.read_value( file_unit, dtype='string' )
        self.BOTTOM         = (self.BC_method.upper() == 'BOTTOM')
        self.RIGHT          = (self.BC_method.upper() == 'RIGHT')
        self.CORNER         = (self.BC_method.upper() == 'CORNER')
        self.FOUR_SIDES     = (self.BC_method.upper() == 'FOUR_SIDES')
        if not(numpy.any([self.BOTTOM, self.RIGHT,
                          self.CORNER, self.FOUR_SIDES])):
            self.BOTTOM = True
        #----------------------------------------------------------------------------
        self.FILL_PITS      = cfg.read_value( file_unit, dtype='boolean' )
        self.ADAPTIVE_DT    = cfg.read_value( file_unit, dtype='boolean' )
        #----------------------------------------------------------------------------
##        print 'BC_method      =', self.BC_method
##        print 'FILL_PITS      =', self.FILL_PITS
##        print 'ADAPTIVE_DT    =', self.ADAPTIVE_DT
        #-----------------------------------------------------------------
        save_grid_dt_info         = cfg.read_input_option( file_unit )
        save_z_grids,  z_gs_file  = cfg.read_output_option( file_unit )
        save_S_grids,  S_gs_file  = cfg.read_output_option( file_unit )
        save_A_grids,  A_gs_file  = cfg.read_output_option( file_unit )
        save_Q_grids,  Q_gs_file  = cfg.read_output_option( file_unit )
        save_Qs_grids, Qs_gs_file = cfg.read_output_option( file_unit )
        #-----------------------------------------------------------------
        self.save_grid_dt  = save_grid_dt_info[1]
        self.SAVE_Z_GRIDS  = save_z_grids
        self.SAVE_S_GRIDS  = save_S_grids
        self.SAVE_A_GRIDS  = save_A_grids
        self.SAVE_Q_GRIDS  = save_Q_grids
        self.SAVE_QS_GRIDS = save_Qs_grids
        self.z_gs_file     = z_gs_file
        self.S_gs_file     = S_gs_file
        self.A_gs_file     = A_gs_file
        self.Q_gs_file     = Q_gs_file
        self.Qs_gs_file    = Qs_gs_file
        #------------------------------------------------------------------
        save_pixels_dt_info        = cfg.read_input_option( file_unit )
        save_z_pixels,  z_ts_file  = cfg.read_output_option( file_unit )
        save_S_pixels,  S_ts_file  = cfg.read_output_option( file_unit )
        save_A_pixels,  A_ts_file  = cfg.read_output_option( file_unit )
        save_Q_pixels,  Q_ts_file  = cfg.read_output_option( file_unit )
        save_Qs_pixels, Qs_ts_file = cfg.read_output_option( file_unit )
        #------------------------------------------------------------------
        self.save_pixels_dt = save_pixels_dt_info[1]
        self.SAVE_Z_PIXELS  = save_z_pixels
        self.SAVE_S_PIXELS  = save_S_pixels
        self.SAVE_A_PIXELS  = save_A_pixels
        self.SAVE_Q_PIXELS  = save_Q_pixels
        self.SAVE_QS_PIXELS = save_Qs_pixels
        self.z_ts_file      = z_ts_file
        self.S_ts_file      = S_ts_file
        self.A_ts_file      = A_ts_file
        self.Q_ts_file      = Q_ts_file
        self.Qs_ts_file     = Qs_ts_file
        
        #-----------------------
        # Close the config file
        #-----------------------
        file_unit.close()

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = maximum(self.save_grid_dt,   self.dt)
        self.save_pixels_dt = maximum(self.save_pixels_dt, self.dt)
        
    #   read_config_file()
    #-------------------------------------------------------------------
    def embed_child_components(self):

        #----------------------------------------------
        # Instantiate and embed "process components"
        # in the place of the CCA ports.
        #----------------------------------------------
        import basins
        import met_base
        ## import uplift_base
        
        self.bp = basins.basins_component()
        self.mp = met_base.met_component()
        ## self.up = uplift_base.uplift_component()
        
    #   embed_child_components()
    #-------------------------------------------------------------------
    def add_child_ports(self):

        if not(self.SILENT):
            print '#########################################'
            print 'NOTE: erosion_base.add_child_ports()'
            print '      is not totally ready yet.'
            print '#########################################'
        self.add_child_port('mp', 'bp')
        
##        self.add_child_port('mp', 'sp')

        ## self.add_child_port('mp', 'pp')   ##### in future ?
        
    #   add_child_ports()
    #-------------------------------------------------------------------
    def initialize_ports(self):

        if not(self.SILENT):
            print '#########################################'
            print 'NOTE: erosion_base.initialize_ports()'
            print '      is not ready yet.'
            print '#########################################'
        return
    
        #-------------------------------------------------
        # Initialize the process objects/components
        # This is also where output files are opened.
        #-------------------------------------------------
        # bp must be initialized first, since other inits
        # use the outlet_ID, etc.
        #--------------------------------------------------
        DEBUG = True
        if (self.bp.get_status() != 'initialized'):        # basin vars
            self.bp.initialize(directory=self.directory,
                               data_prefix=self.data_prefix,
                               case_prefix=self.case_prefix)
            if (DEBUG): print '\nEROSION component initialized BASINS.'
            
        if (self.mp.get_status() != 'initialized'):      # met vars + precip 
            self.mp.initialize(directory=self.directory,
                               data_prefix=self.data_prefix,
                               case_prefix=self.case_prefix)
            if (DEBUG): print '\nEROSION component initialized METEOROLOGY.'

        #-----------------
        # For future use
        #-----------------
##        if (self.up.get_status() != 'initialized'):      # uplift_vars 
##            self.up.initialize(directory=self.directory,
##                               data_prefix=self.data_prefix,
##                               case_prefix=self.case_prefix)
##            if (DEBUG): print '\nEROSION component initialized UPLIFT.'
            
    #   initialize_ports()
    #-------------------------------------------------------------------
    def initialize_d8_vars(self):

        #---------------------------------------------
        # Compute and store a variety of (static) D8
        # flow grid variables.  Embed structure into
        # the "channel_base" component.
        #---------------------------------------------
        self.d8 = d8_base.d8_component()
        self.d8.initialize(directory=self.directory,
                           data_prefix=self.data_prefix,
                           case_prefix=self.case_prefix,
                           SILENT=self.SILENT,
                           REPORT=self.REPORT)

        #--------------------------------------        
        # This is called by update() function
        #--------------------------------------
##        self.d8.update(self.time, SILENT=False, REPORT=True)

    #   initialize_d8_vars()
    #-------------------------------------------------------------
    def initialize_DEM(self, SILENT=True):

        #-------------------------------------------------------
        # Notes: This function initializes the DEM either from
        #        self.z0_file or using create_initial_DEM().
        #-------------------------------------------------------
        ### if (self.z0_file != ''):
        if (self.make_z0_method == 'READ_FILE'):
            #---------------------------------------
            # Read inital elevation grid from file
            #---------------------------------------
            DEM_unit  = open(self.z0_file, 'rb')
            file_size = os.path.getsize(DEM_unit.name)
            dbl_size  = self.nx * self.ny * int32(8)
            flt_size  = self.nx * self.ny * int32(4)
            int_size  = self.nx * self.ny * int32(2)
            if (file_size == dbl_size):
                RTG_type = 'DOUBLE'
            elif (file_size == flt_size):
                RTG_type = 'FLOAT'
            elif (file_size == int_size):
                RTG_type = 'INTEGER'
            else:
                print 'ERROR in initialize_DEM().'
                print '   Cannot determine DEM data type.'
                return
            
            self.DEM = rtg_files.read_grid( self.z0_file, self.rti,
                                            RTG_type=RTG_type ) 
            if (RTG_type != 'FLOAT'):
                self.DEM = float32(self.DEM)
        else:    
            #--------------------------------
            # Create initial elevation grid
            #--------------------------------
            self.create_initial_DEM()

            #-------------------------------------
            # See self.update_DEM_edge_values().
            #----------------------------------------------------------
            # 01/13/06. Otherwise we end up with very large values on
            # the indicated boundary that prevent good use of color
            #----------------------------------------------------------
            # self.update_DEM_edge_values()

            #------------------------------------
            # Save the initial DEM to a file ??
            #------------------------------------
            #### Check for overwrite here !!!  ##############################
##            DEM_unit = open(self.z0_file, 'wb')
##            if (self.rti.SWAP_ENDIAN):   ## (rti undefined so far)
##                self.DEM.byteswap(True)
##            self.DEM.tofile(DEM_unit)
##            DEM_unit.close()

        zmin = numpy.nanmin( self.DEM )
        zmax = numpy.nanmax( self.DEM )
        if not(SILENT):
            print 'Initial (z_min, z_max) =', zmin, zmax
            
    #   initialize_DEM()
    #-------------------------------------------------------------
    def create_initial_DEM(self):

        #------------------------------------------------------
        # Notes: This routine allows a "z0_method" and a
        #        "noise_method" to be combined to generate
        #        an initial surface (z0) with optional noise.
        #
        #        z0_methods:    FLAT, PLANE, CORNER_PLANE
        #        noise_methods: GAUSSIAN, MIDPOINTS

        #        Noise grid is generated first, which is
        #        more efficient in the FLAT case.

        #        "seed" is read from CFG file and should be
        #        a 4- or 5-digit integer (e.g. 36421).
        #        If the same seed is used, the same sequence
        #        of random numbers is generated, which allows
        #        for reproducible results and comparisons.
        #------------------------------------------------------
        nx = self.nx  # (local synonyms)
        ny = self.ny

        #---------------------------------------
        # Uncorrelated Gaussian random numbers
        # which essentially models white noise
        #---------------------------------------
        # When sigma or factor = 1, then range
        # is pretty much -3 to 3.
        #---------------------------------------
        if (self.GAUSSIAN):    
            numpy.random.seed( self.seed )
            self.DEM = random.normal(loc=0.0, scale=1.0, size=(ny, nx))
               #(mean = 0.0, stddev = 1.0)
            #-----------------------------------
            # factor = (1 / float64(3))               #(-3,3) -> (-1,1) 
            # factor = factor * float64(300)
            factor = self.noise_scale
            
            #-----------------------------------
            # Slope should dominate over noise
            #-----------------------------------
            #*** factor = factor * (slope / 2d)
            self.DEM = factor * self.DEM
            #------------------------------
            # if (PLANE OR CORNER_PLANE) then mean=0.0 else mean=2.0
            # stddev = 4.0  ;(1.0)
            # DEM = (stddev * DEM) + mean
        elif (self.MIDPOINTS):
            #-----------------------------------------------------------
            # Use "midpoint displacement" to create a fractal surface
            # with correlation structure as in MARSSIM (A. Howard).
            # See midpoints.py in code directory.
            #-----------------------------------------------------------
            nn = max(nx, ny)
            n_levels = numpy.ceil(numpy.log(nn-1) / numpy.log(2))
            n_levels = int16( n_levels )
            surf = midpoints.make_fractal_surface( n_levels, H=1.5,
                                                   scale=self.noise_scale,
                                                   seed=self.seed,
                                                   SILENT=True)
            self.DEM = surf[0:ny, 0:nx]
        else:    
            self.DEM = float32(0)
         
        #--------------------------------------
        # Inclined plane tilted toward bottom
        #--------------------------------------
        if (self.PLANE):
            ## IDs = reshape(arange(nx*ny, dtype='Int32'), [ny, nx])
            IDs  = self.d8.ID_grid
            cols = (IDs % nx)
            rows = (IDS / nx)
            x    = (self.dx * cols)
            y    = (self.dy * (ny - rows))
            #-------------------------------------
            z    = (self.z0_plane_dz_dx * x) + \
                   (self.z0_plane_dz_dy * y)
            #-------------------------------------
            self.DEM += z
        
        #-------------------------------------------------
        # Inclined plane tilted toward lower left corner
        #-------------------------------------------------
        if (self.CORNER_PLANE):
            ## IDs = reshape(arange(nx*ny, dtype='Int32'), [ny, nx])
            IDs  = self.d8.ID_grid
            cols = (IDs % nx)
            rows = (IDS / nx)
            x    = (self.dx * cols)
            y    = (self.dy * (ny - rows))
            #--------------------------------------------
            a = (self.z0_plane_S / sqrt(float32(2.0)))
            b = a
            z = (a * x) + (b * y)
            #--------------------------------------------
##            z    = (self.z0_plane_dz_dx * x) + \
##                   (self.z0_plane_dz_dy * y)
            #-------------------------------------
            self.DEM += z

        #----------------------------
        # Make sure type is FLOAT ?
        #----------------------------
        self.DEM = numpy.float32(self.DEM)
        
        #-------------------------
        # Save new DEM to a file
        #------------------------- 
        rtg_files.write_grid( self.DEM, self.z0_file, self.rti)

    #   create_initial_DEM()
    #-------------------------------------------------------------
    def initialize_boundary_conditions(self):

        nx = self.nx   # (local synonyms)
        ny = self.ny
        ID_type = 'Int32'
        
        #------------------------------------
        # Use bottom row/edge as base level
        #------------------------------------
        if (self.BOTTOM):    
            self.base_IDs = arange(nx, dtype=ID_type) + (nx * (ny - 1))
            #*** above_row = (base_IDs - nx)
            #*** base_IDs  = [above_row, base_IDs]
            
            #---------------------------------
            # Change values in the top row ?
            #----------------------------------
            # Copy values so slope & fluxes
            # are zero at top edge for PLANE1
            #----------------------------------
            top_IDs  = arange(nx, dtype=ID_type)
            row2_IDs = arange(nx, dtype=ID_type) + nx
            self.DEM[top_IDs] = self.DEM[row2_IDs]
            #*** self.DEM[top_IDs] = 0.0
        
        #-------------------------------
        # Use right edge as base level
        #-------------------------------
        if (self.RIGHT):    
            self.base_IDs = nx * (arange(ny, dtype=ID_type) + 1)
            self.base_IDs = self.base_IDs - 1
            #*** prev_col = (base_IDs - 1L)
            #*** base_IDs = [prev_col, base_IDs]
            
            #---------------------------------
            # Change values on the left edge
            #---------------------------------
            left_IDs = arange(ny, dtype=ID_type) * nx
            #** emax = max(DEM, /NAN)
            #** DEM[left_IDs] = (emax * 0.2)
            self.DEM[left_IDs] = float32(0.0)
        
        #-----------------------------------
        # Use all four sides as base level
        #-----------------------------------
        if (self.FOUR_SIDES):    
            T_IDs = arange(nx, dtype=ID_type)
            B_IDs = T_IDs + nx * (ny - 1)
            L_IDs = nx * (arange(ny - 2, dtype=ID_type) + 1)
            R_IDs = L_IDs + (nx - 1)
            self.base_IDs = concatenate((T_IDs, B_IDs, L_IDs, R_IDs))
        
        #---------------------------------------
        # Use bottom left corner as base level
        #---------------------------------------
        if (self.CORNER):    
            ID1 = nx * (ny - 1)    # (lower left pixel)
            ID2 = ID1 - nx         # (just above ID1)
            ID3 = ID1 - (2 * nx)   # (just above ID2)
            self.base_IDs = concatenate(( ID1, ID1 + 1, ID1 + 2,
                                          ID2, ID2 + 1, ID2 + 2,
                                          ID3, ID3 + 1, ID3 + 2 ))
        
        #--------------------------------------
        # Set the initial base-level height ?
        #--------------------------------------
        if (self.BOTTOM or self.RIGHT or self.CORNER or self.FOUR_SIDES):    
            #-------------------------------------
            # Subtracting 1 here is not good for
            # continuing on from a previous run
            #-------------------------------------
            self.base_level = nanmin(self.DEM)
        else:    
            self.base_level = float32(0)
        
    #   initialize_boundary_conditions()  
    #-------------------------------------------------------------------
    def initialize_computed_vars(self):

        self.dt    = float32(0.1)  # [years]
        
        self.vol_R = float64(0)  # (for mass balance)
        self.vol_U = float64(0)  # (for mass balance)
        
        self.dz_max_vec = zeros([self.n_steps], dtype='Float32')
        self.dz_max     = float32(-9999)
        
    #   initialize_computed_vars()
    #-------------------------------------------------------------------
    def update_R(self):

        #------------------------------------------------------
        # Note: self.R is currently set by read_config_file()
        #------------------------------------------------------
        return

        ##################################################
        ##################################################
        ##  CONVERT UNITS FROM [m/s] to [m/yr] BELOW !!
        ##################################################
        ##################################################
    
        #----------------------------------------
        # Compute the "excess rainrate", R.
        # Each term must have same units: [m/s]
        # Sum = net gain/loss rate over pixel.
        #----------------------------------------------------
        # R can be positive or negative.  If negative, then
        # water is removed from the surface at rate R until
        # surface water is consumed.
        #--------------------------------------------------------------
        # P  = precip_rate   [m/s]  (converted by read_input_data()).
        # SM = snowmelt rate [m/s]
        # GW = seep rate     [m/s]  (water_table intersects surface)
        # ET = evap rate     [m/s]
        # IN = infil rate    [m/s]
        # MR = icemelt rate  [m/s]
        #--------------------------------------------------------------        
        P  = self.get_port_data('P',  self.mp,  'METEOROLOGY')
        # SM = self.get_port_data('SM', self.sp,  'SNOW')
        # GW = self.get_port_data('GW', self.gp,  'SATZONE')
        # ET = self.get_port_data('ET', self.ep,  'EVAP')
        # IN = self.get_port_data('IN', self.ip,  'INFIL')
        # MR = self.get_port_data('MR', self.iip, 'ICE')
        
        #--------------
        # For testing
        #--------------        
##        print '(Pmin,  Pmax)  =', P.min(),  P.max()
##        print '(SMmin, SMmax) =', SM.min(), SM.max()
##        print '(GWmin, GWmax) =', GW.min(), GW.max()
##        print '(ETmin, ETmax) =', ET.min(), ET.max()
##        print '(INmin, INmax) =', IN.min(), IN.max()
##        print '(MRmin, MRmax) =', MR.min(), MR.max()
##        # print '(Hmin,  Hmax)  =', H.min(), H.max()
##        print ' '

        self.R = P
        ## self.R = (P + SM + GW + MR) - (ET + IN)
            
    #   update_R()
    #-------------------------------------------------------------------
    def update_R_integral(self):

        #-----------------------------------------------
        # Update mass total for R, sum over all pixels
        #-----------------------------------------------   
        volume = double(self.R * self.da * self.dt)  # [m^3]
        if (size(volume) == 1):
            self.vol_R += (volume * self.rti.n_pixels)
        else:
            self.vol_R += sum(volume)

    #   update_R_integral()
    #-------------------------------------------------------------------
    def update_U(self):

        #------------------------------------------------------
        # Note: self.U is currently set by read_config_file()
        #------------------------------------------------------
        return

        #----------------------------------------
        # Compute the "uplift rate", U, [mm/yr]
        #----------------------------------------      
        self.U = self.get_port_data('U',  self.up,  'UPLIFT')
        
    #   update_U()     
    #-------------------------------------------------------------------
    def update_U_integral(self):

        #-----------------------------------------------
        # Update mass total for U, sum over all pixels
        #-----------------------------------------------   
        volume = double(self.U * self.da * self.dt)  # [m^3]
        if (size(volume) == 1):
            self.vol_U += (volume * self.rti.n_pixels)
        else:
            self.vol_U += sum(volume)

    #   update_U_integral() 
    #-------------------------------------------------------------------
    def update_base_level(self):
        
        #--------------------------------------
        # Lower base level for bottom row or
        # rightmost column or LL corner, etc.
        #--------------------------------------
        if (self.BOTTOM or self.RIGHT or self.CORNER or self.FOUR_SIDES):    
            #---------------------------------------
            # NB!  Inside loop since dt is dynamic
            #---------------------------------------
            # Units of BLR are mm/year so we must
            # convert mm/yr to meters/yr.
            #---------------------------------------
            drop = (self.BLR / self.mm_per_m * self.dt)
            self.base_level -= drop
            self.DEM.flat[ self.base_IDs ] = self.base_level
            
        #-------------------------------
        # Maintain boundary condition   (Used for Test5)
        # e.g. zero out all four edges
        #-------------------------------
        # nx = self.nx
        # ny = self.ny
        # self.DEM[0,:]     = 0.0    # (left edge)
        # self.DEM[nx-1, :] = 0.0    # (right edge)
        # self.DEM[:, 0]    = 0.0    # (top edge)
        # self.DEM[:, ny-1] = 0.0    # (bottom edge)
    
    #   update_base_level()
    #-------------------------------------------------------------------
    def update_DEM_edge_values(self):
        
        #-------------------------------------------
        # 01/16/06.  Adjust DEM edge values since
        # they can't erode and will otherwise stay
        # big and skew the color stretches.
        #-------------------------------------------
        if (self.BOTTOM):
            self.DEM[0,:] = self.DEM[1,:] + float32(0.01)
        if (self.RIGHT):
            self.DEM[:,0] = self.DEM[:,1] + float32(0.01)

##        if (self.BOTTOM):
##            self.DEM[:, 1] = numpy.nanmax( self.DEM )
##            self.DEM[:, 0 ]= self.DEM[:, 1] + 0.01
##        if (self.RIGHT):
##            self.DEM[1, :] = numpy.nanmax( self.DEM )
##            self.DEM[0, :] = self.DEM[1, :] + 0.01

    #   update_DEM_edge_values()
    #-------------------------------------------------------------------
    def update_d8_vars(self, SILENT=True, REPORT=False):

        #---------------------------------------------
        # Update the D8 flow grid and all vars that
        # depend on it, including D8 area grid.
        #---------------------------------------------
        # Area grid units are either 'm^2' or 'km^2'
        # based on a setting in "*_d8.cfg" file.
        # All length units are given in meters.
        #---------------------------------------------
        # d8.update() needs a depression-filled DEM
        # and can later get it from a CCA port.
        #---------------------------------------------        
        self.d8.update( self.time, self.DEM,
                        SILENT=SILENT, REPORT=REPORT)
        
    #   update_d8_vars()
    #-------------------------------------------------------------------
    def update_slope_grid(self, SILENT=True, REPORT=False):

        #----------------------------------------------------
        # Notes: Make sure that d8 component is initialized
        #        with the same directory, data_prefix, etc.
        #        as this component.  Otherwise, the "shape"
        #        of DEM and ds, A, etc. won't match.
        #----------------------------------------------------
        if not(SILENT):    
            print 'Updating slope grid...'
   
        #--------------------------------------------------
        # Compute slope (rise/run) toward D8 parent pixel
        #--------------------------------------------------
        # pIDs gives indices in the "numpy.where" style
        #--------------------------------------------------        
        pIDs   = self.d8.parent_IDs
        self.S = ((self.DEM - self.DEM[pIDs]) / self.d8.ds)

##        w = where(self.DEM < self.DEM[pIDs])
##        n_bad = size(w[0])
##        if (n_bad != 0):
##            print '   Number of uphill parents =', n_bad
            
        #---------------------------------------------------
        # Set slope to zero wherever the D8 flow direction
        # is undefined.  This makes Qs(Q,S) = 0.
        #---------------------------------------------------
        ## w  = where(self.d8.parent_ID_grid == 0)
        w  = where(self.d8.flow_grid == 0)
        nw = size(w[0])
        if (nw != 0):    
            self.S[w] = 0

        #-------------------------
        # Check for bad S values
        #-------------------------
        w2  = where(self.S < 0)
        nw2 = size(w2[0])
        if (nw2 != 0):    
            # self.S[w2] = 1e-8
            print '   Negative values found in slope grid.'
            print '   Found ' + str(nw2) + ' negative values.'
            
            #---------------------------
            # Return the cols and rows
            #---------------------------
            rows = w2[0]
            cols = w2[1]
            print '   cols ='
            print cols
            print '   rows ='
            print rows
            sys.exit()   #####################
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            ## S_str = str(self.S.min())  + ', '  + str(self.S.max())
            S_str = str(nanmin(self.S)) + ', ' + str(nanmax(self.S))
            print '    min(S), max(S) = ' + S_str + ' [m/m]'

    #   update_slope_grid()
    #-------------------------------------------------------------------
    def update_Q_grid(self, SILENT=True, REPORT=False):

        #--------------------------------------------------------------
         #NOTES:  Q = annual discharge [m^3 / yr]
        #         R = geomorphically effective rainrate [meters / yr]
        #             (unless p ne 0, then R = coefficient)
        #         A = contributing area [meters^2]
        #--------------------------------------------------------------
        if not(SILENT):    
            print 'Updating discharge grid...'
        
        if (self.p != 1):    
            self.Q = self.R * (self.d8.A ** self.p)
        else:    
            #----------------------
            # Don't convert units
            #----------------------
            #** R2 = R / self.secs_per_year   ;[m/yr -> m/s]
            self.Q = self.R * self.d8.A
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            ## Q_str = str(self.Q.min())  + ', '  + str(self.Q.max())
            Q_str = str(nanmin(self.Q)) + ', ' + str(nanmax(self.Q))
            print '    min(Q), max(Q) = ' + Q_str + ' [m^3/yr]'

    #   update_Q_grid()
    #-------------------------------------------------------------------
    def update_Qs_grid(self, SILENT=True, REPORT=False):

        #--------------------------------------------------------
        #NOTES:  Qs = annual sed. discharge [m^3/ yr]
        #        Q  = annual discharge [m^3 / yr]
        #        S  = slope [unitless]
        #        k  = coefficient [(m^3 / yr)^(1 - m)]
        #        m  = discharge exponent (usually in [1,2])
        #        n  = slope exponent (usually in [1,2])

        #        The standard formula Qs = k * Q^m * S^n is for
        #        the case where Q and Qs have units of m^3/sec.
        #        The extra factor below makes the formula valid
        #        for the case where units of both are m^3/yr.
        #        That is, Qs' = fac * kf * (Q')^m * S^n.

        #        The Slope_Grid function checks for negative
        #        slopes and either adjusts them or aborts.
        #--------------------------------------------------------
        if not(SILENT):
            print 'Updating sed. discharge grid...'
        
        fac1    = self.secs_per_year ** (float64(1) - self.m)
        fac2    = self.K * (self.Q ** self.m)
        fac3    = (self.S ** self.n)
        self.Qs = fac1 * fac2 * fac3     #[m^3 / year]
        
        #---------------------------------------
        # Impose no-flux boundary condition
        # on all four edges of DEM ?
        # NB!  If S=0 on edges, then Qs=0.
        #---------------------------------------
        # Is it better to put a "wall" around
        # the outside to impose this B.C. ?
        #---------------------------------------
        # nx = self.nx
        # ny = self.ny
        # self.Qs[0, :]    = 0.0
        # self.Qs[nx-1, :] = 0.0
        # self.Qs[:, 0]    = 0.0
        ### self.Qs[:, ny-1] = 0.0    ;(Let sediment get out)
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):    
            ## Qs_str = str(self.Qs.min())  + ', '  + str(self.Qs.max())
            Qs_str = str(nanmin(self.Qs)) + ', ' + str(nanmax(self.Qs))
            print '    min(Qs), max(Qs) = ' + Qs_str + ' [m^3/yr]'

    #   update_Qs_grid()
    #-------------------------------------------------------------------
    def update_DEM(self, SILENT=True, REPORT=False):

        #------------------------------------------------------------
        # Notes: This version computes dt dynamically to ensure
        #        stability.

        #        DEM   = current elevation grid  [m]
        #        Qs    = kf * Q^mf * S^nf = sed. discharge [m^3/yr]
        #        da    = pixel area grid  [m^2]
        #        dt    = channel flow timestep  [s]  (returned)
        #        del_z = elevation drops to parent pixels [m]
        #        dz    = net elevation change due to sed. flux [m]
        #        U     = tectonic uplift rate [mm/year]
        #        w1    = IDs of pixels that...
        #        p1    = IDs of parent pixels that...

        #        NB!  Don't want elevations to change at pixels
        #        where base level is fixed, as with flow to sea.
        #------------------------------------------------------------
        if not(SILENT):    
            print 'Updating elevations...'
        
        #-----------------------------
        # Initialize dQs with Qs_out
        #-----------------------------
        dQs = -self.Qs + (self.U * (self.da / self.mm_per_m))   #(m^3 per year)
        
        #---------------------------------------------
        # Add contributions from neighbor pixels
        # This approach assumes fixed-length pixels.
        #---------------------------------------------
        if (self.d8.n1 != 0):    
            dQs[self.d8.p1] += self.Qs[self.d8.w1]
        if (self.d8.n2 != 0):    
            dQs[self.d8.p2] += self.Qs[self.d8.w2]
        if (self.d8.n3 != 0):    
            dQs[self.d8.p3] += self.Qs[self.d8.w3]
        if (self.d8.n4 != 0):    
            dQs[self.d8.p4] += self.Qs[self.d8.w4]
        if (self.d8.n5 != 0):    
            dQs[self.d8.p5] += self.Qs[self.d8.w5]
        if (self.d8.n6 != 0):    
            dQs[self.d8.p6] += self.Qs[self.d8.w6]
        if (self.d8.n7 != 0):    
            dQs[self.d8.p7] += self.Qs[self.d8.w7]
        if (self.d8.n8 != 0):    
            dQs[self.d8.p8] += self.Qs[self.d8.w8]
        
        #---------------------------
        # Compute downstream drops
        #---------------------------
        pIDs  = self.d8.parent_IDs
        del_z = numpy.maximum((self.DEM - self.DEM[ pIDs ]), 0)
        
        #------------------------
        # Experiment 2: 8/10/04
        #------------------------
        del_Qs = numpy.maximum((self.Qs - self.Qs[pIDs]), 0)
        
        #------------------------
        # Experiment 1: 8/10/04
        #------------------------
        # nx  = self.nx
        # ny  = self.ny
        # dgm = zeros((ny,nx), dtype='Float32') + 99999.0
        #---------------------------------------------------
        # if (n1 != 0): dgm[p1] = minimum(dgm[p1], dg[w1])
        # if (n2 != 0): dgm[p2] = minimum(dgm[p2], dg[w2])
        # if (n3 != 0): dgm[p3] = minimum(dgm[p3], dg[w3])
        # if (n4 != 0): dgm[p4] = minimum(dgm[p4], dg[w4])
        # if (n5 != 0): dgm[p5] = minimum(dgm[p5], dg[w5])
        # if (n6 != 0): dgm[p6] = minimum(dgm[p6], dg[w6])
        # if (n7 != 0): dgm[p7] = minimum(dgm[p7], dg[w7])
        # if (n8 != 0): dgm[p8] = minimum(dgm[p8], dg[w8])
        #---------------------------------------------------
        # wn  = where(dgm >= 99999.0)
        # nwn = size(wn[0])
        # if (nwn != 0): dgm[wn]=0.0  ;*******************
        # dg = dgm  ;*********************
        
        #----------------------------------
        # Compute largest stable timestep
        # with a small factor of safety
        #---------------------------------------------
        # Decreased factor from 0.8 to 0.4 (8/10/04)
        #---------------------------------------------
        # Tried factor of 0.6 on 8/12/04, but for
        # mf = nf = 1, one or more pits would form
        # near the outlet so decreased back to 0.4.
        #---------------------------------------------
        # Set dt so that dz < dg
        #--------------------------------
        # del_Qs = dQs  ;**********
        wp  = where(logical_and((del_Qs > 0), (del_z > 0)))
        nwp = size(wp[0])    # (deposition sites)
        if (nwp != 0):
            if (size(self.da) == 1):
                vals = self.da * (del_z[wp] / del_Qs[wp])
            else:
                vals = self.da[wp] * (del_z[wp] / del_Qs[wp])               
            self.dt = float64(0.4) * numpy.nanmin( vals )  # [years]
        else:
            self.dt = float64(1) # [years]
        
        #-----------------------------------
        # Should there be a max timestep ?
        #-----------------------------------
        #*** dt_max = 5d
        #*** self.dt = (self.dt < dt_max)
        #*** print,'dt = ' + str(self.dt)
        
        #-------------------------------
        # Don't let dt get too small ?
        #-----------------------------------
        # dt_min must be less than 1e-4
        # for case mf=1.5, nf=1.0, kf=1.0.
        # Making kf smaller allows dt_min
        # to be bigger. Now kf is smaller.
        #-----------------------------------
        # Fixed units issue in Qs_Grid so
        # should be able to reduce dt_min.
        #-----------------------------------
        dt_min = float32(1E-2)
        if (self.dt < dt_min):    
            print '******************************************'
            print 'Aborting: Stable dt is too small.'
            print 'Computed dt = ' + str(self.dt)
            print '******************************************'
            print ' '
            sys.exit()
        
        #-----------------------------------------------------------
        # Compute dz and the largest dz obtained during this
        # timestep. This is used to check stability & convergence.
        #-----------------------------------------------------------
        dz = dQs * (self.dt / self.da)
        self.dz_max = nanmax(dz)
        self.dz_max_vec[ self.time_index - 1 ] = self.dz_max
        
        #------------------
        # Optional report
        #------------------
        if (REPORT):
            dz_str     = str(nanmin(dz))   + ', ' + str(self.dz_max)
            del_z_str  = str(del_z.min())  + ', ' + str(del_z.max())
            del_Qs_str = str(del_Qs.min()) + ', ' + str(del_Qs.max())
            print '    min(dz),     max(dz)     = ' + dz_str + ' [m]'
            print '    min(del_z),  max(del_z)  = ' + del_z_str
            print '    min(del_Qs), max(del_Qs) = ' + del_Qs_str
        
        #----------------
        # Add dz to DEM
        #----------------
        self.DEM += float32(dz)
        self.dz   = float32(dz)   # (save for retrieval by caller ?)
 
    #   update_DEM()
    #-------------------------------------------------------------------
    def fill_pits_in_DEM(self, SILENT=True):

        if not(SILENT):    
            print 'Filling depressions in DEM...'

##        DEM_before = self.DEM.copy()  # (For testing)
        
        fill_pits.fill_pits(self.DEM, 'FLOAT', self.nx, self.ny,
                            SILENT=SILENT)

##        #--------------
##        # For testing
##        #--------------
##        w  = where(DEM_before != self.DEM)
##        nw = size(w[0])
##        print 'Number of pixels changed by fill_pits =', nw
        
    #   fill_pits_in_DEM()
    #-------------------------------------------------------------------
    def update_DEM_min_and_max(self, REPORT=False):
        
        self.DEM_min = numpy.nanmin( self.DEM )
        self.DEM_max = numpy.nanmax( self.DEM )

        #------------------
        # Optional report
        #------------------
        if (REPORT):
            z_str = str(self.DEM_min) + ', ' + str(self.DEM_max)
            print '    min(z), max(z) = ' + z_str + ' [m]'
            
    #   update_DEM_min_and_max() 
    #-------------------------------------------------------------------
    def check_stability(self):
        
        #----------------------------------
        # Check for one type of stability
        # (This may be obsolete now.)
        #----------------------------------
        if (self.dz_max > float32(200)):    
            print '************************************************'
            print 'Program aborted due to very large change.'
            print 'Time step or K is probably too large'
            print 'for this set of input parameters.'
            print '   dx = ' + str(self.dx) + ' [meters]'
            print '   dy = ' + str(self.dy) + ' [meters]'
            print '   dt = ' + str(self.dt) + ' [years]'
            print '   K  = ' + str(self.K)
            print '   m  = ' + str(self.m)
            print '   n  = ' + str(self.n)
            print '************************************************'
            print ' '
            return False
            ## sys.exit()
        else:
            return True
        
    #   check_stability()
    #-------------------------------------------------------------------
    def check_finished(self):

        #---------------------------------------------------------
        # Note: TINY_DZ can occur either because dt required for
        #       stability is really small or because we have
        #       converged to a steady-state landscape.
        #---------------------------------------------------------
        #       If self.DONE has already been set to True by
        #       another function or component, this function
        #       preserves that setting (see below).
        #---------------------------------------------------------        
        TIMES_UP  = (self.time_index >= self.n_steps)
        TINY_DZ   = (self.dz_max < self.dz_tolerance)
        self.DONE = (self.DONE or TIMES_UP or TINY_DZ)
        
        if (TINY_DZ):
            tol_str = str(self.dz_tolerance)
            print 'Aborting since dz_max < ' + tol_str + '.'
                     
    #   check_finished()
    #-------------------------------------------------------------------
##    def check_steady_state(self):
##
##    #   check_steady_state()
    #-------------------------------------------------------------------
    def open_input_files(self):

        if (self.make_z0_method == 'READ_FILE'):
            self.z0_unit = model_input.open_file(self.z0_type, self.z0_file)

    #   open_input_files()        
    #-------------------------------------------------------------------  
    def read_input_files(self):

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        if (self.make_z0_method == 'READ_FILE'):
            self.z0 = model_input.read_next(self.z0_unit, self.z0_type, self.rti)
            self.z0_unit.close()   #########
        
##        slopes = model_input.read_next(self.slope_unit, self.slope_type, rti)
##        if (slopes != None): self.slopes = slopes

    #   read_input_files()     
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.make_z0_method == 'READ_FILE'):
            self.z0_unit.close()
            ### if (self.z0_file != ''): self.z0_unit.close()

    #   close_input_files()       
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        #--------------------------------------------------------
        # Notes:  Whenever the case_prefix changes (e.g. user
        #         changes it and clicks the Next button in the
        #         "run var" panel), all output filenames should
        #         be reset to new defaults.
        #--------------------------------------------------------
        self.z_gs_file  = (self.case_prefix + '_2D-z.rts')
        self.S_gs_file  = (self.case_prefix + '_2D-S.rts')
        self.A_gs_file  = (self.case_prefix + '_2D-A.rts')
        self.Q_gs_file  = (self.case_prefix + '_2D-Q.rts')
        self.Qs_gs_file = (self.case_prefix + '_2D-Qs.rts')
        #----------------------------------------------------
        self.z_ts_file  = (self.case_prefix + '_0D-z.txt')
        self.S_ts_file  = (self.case_prefix + '_0D-S.txt')
        self.A_ts_file  = (self.case_prefix + '_0D-A.txt')
        self.Q_ts_file  = (self.case_prefix + '_0D-Q.txt')
        self.Qs_ts_file = (self.case_prefix + '_0D-Qs.txt')
        
    #   update_outfile_names()     
    #-------------------------------------------------------------------  
    def open_output_files(self):
         
        #--------------------------------------
        # Open new files to write grid stacks
        #--------------------------------------
        if (self.SAVE_Z_GRIDS):   
            self.z_gs_unit = open_new_gs_file(self.z_gs_file, self.rti,
                                              grid_name='z',
                                              long_name='elevation grid',
                                              units_name='m')
            
        if (self.SAVE_S_GRIDS):    
            self.S_gs_unit = open_new_gs_file(self.S_gs_file, self.rti,
                                              grid_name='S',
                                              long_name='slope_grid',
                                              units_name='m/m')
        
        if (self.SAVE_A_GRIDS):    
            self.A_gs_unit = open_new_gs_file(self.A_gs_file, self.rti,
                                              grid_name='A',
                                              long_name='contributing_area_grid',
                                              units_name='km^2')

        if (self.SAVE_Q_GRIDS):    
            self.Q_gs_unit = open_new_gs_file(self.Q_gs_file, self.rti,
                                              grid_name='Q',
                                              long_name='water_discharge_grid',
                                              units_name='m^3/s')

        if (self.SAVE_QS_GRIDS):    
            self.Qs_gs_unit = open_new_gs_file(self.Qs_gs_file, self.rti,
                                               grid_name='Qs',
                                               long_name='sediment_discharge_grid',
                                               units_name='m^3/s')  ###########
            
        #---------------------------------------
        # Open text files to write time series
        #---------------------------------------
        IDs = self.outlet_IDs
        if (self.SAVE_Z_PIXELS):    
            self.z_ts_unit = open_new_ts_file( self.z_ts_file, 'z', IDs, 'years')
            
        if (self.SAVE_S_PIXELS):
            self.S_ts_unit = open_new_ts_file( self.S_ts_file, 'S', IDs, 'years')
            
        if (self.SAVE_A_PIXELS):    
            self.A_ts_unit = open_new_ts_file( self.A_ts_file, 'A', IDs, 'years')
            
        if (self.SAVE_Q_PIXELS):    
            self.Q_ts_unit = open_new_ts_file( self.Q_ts_file, 'Q', IDs, 'years')

        if (self.SAVE_QS_PIXELS):    
            self.Qs_ts_unit = open_new_ts_file( self.Qs_ts_file, 'Qs', IDs, 'years')
            
    #   open_output_files()
    #-------------------------------------------------------------------  
    def write_output_files(self, time=None):

        #---------------------------------------------------------
        # Notes:  This function was written to use only model
        #         time (maybe from a caller) in seconds, and
        #         the save_grid_dt and save_pixels_dt parameters
        #         read by read_config_file().
        #
        #         read_config_file() makes sure that all of
        #         the "save_dts" are larger than or equal to the
        #         process dt.
        #---------------------------------------------------------
        
        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time is None):
            time = self.time
        model_time = int(time)
        
        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
##        if ((self.time_index % self.grid_save_step) == 0):
##             self.save_grids()
##        if ((self.time_index % self.pixel_save_step) == 0):
##             self.save_pixel_values()
        
    #   write_output_files()
    #-------------------------------------------------------------------  
    def close_output_files(self):

        if (self.SAVE_Z_GRIDS):   self.z_gs_unit.close()   
        if (self.SAVE_S_GRIDS):   self.S_gs_unit.close()    
        if (self.SAVE_A_GRIDS):   self.A_gs_unit.close()   
        if (self.SAVE_Q_GRIDS):   self.Q_gs_unit.close()
        if (self.SAVE_QS_GRIDS):  self.Qs_gs_unit.close()
        #---------------------------------------------------
        if (self.SAVE_Z_PIXELS):  self.z_ts_unit.close()   
        if (self.SAVE_S_PIXELS):  self.S_ts_unit.close()    
        if (self.SAVE_A_PIXELS):  self.A_ts_unit.close()    
        if (self.SAVE_Q_PIXELS):  self.Q_ts_unit.close()
        if (self.SAVE_QS_PIXELS): self.Qs_ts_unit.close()
        
    #   close_output_files()              
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        #---------------------------------------------------
        # Notes:  Each variable is saved as a grid whether
        #         it is a scalar or already a (2D) grid.
        #---------------------------------------------------
        nx = self.nx
        ny = self.ny

        if (self.SAVE_Z_GRIDS):
            save_as_grid_to_file(self.z_gs_unit, self.DEM, 'z', nx, ny)
            
        if (self.SAVE_S_GRIDS):
            save_as_grid_to_file(self.S_gs_unit, self.S, 'S', nx, ny)
            
        if (self.SAVE_A_GRIDS):
            save_as_grid_to_file(self.A_gs_unit, self.d8.A, 'A', nx, ny)

        if (self.SAVE_Q_GRIDS):
            save_as_grid_to_file(self.Q_gs_unit, self.Q, 'Q', nx, ny)    

        if (self.SAVE_QS_GRIDS):
            save_as_grid_to_file(self.Qs_gs_unit, self.Qs, 'Qs', nx, ny)
            
    #   save_grids()
    #-------------------------------------------------------------------  
    def save_pixel_values(self):   ##### save_time_series_data(self)  #######
        
        IDs      = self.outlet_IDs
        
        if (self.SAVE_Z_PIXELS):
            write_ts_file_line(self.z_ts_unit, self.time, self.DEM, IDs)
                    
        if (self.SAVE_S_PIXELS):
            write_ts_file_line(self.S_ts_unit, self.time, self.S, IDs)
            
        if (self.SAVE_A_PIXELS):
            write_ts_file_line(self.A_ts_unit, self.time, self.d8.A, IDs)
            
        if (self.SAVE_Q_PIXELS):
            write_ts_file_line(self.Q_ts_unit, self.time, self.Q, IDs)
            
        if (self.SAVE_QS_PIXELS):
            write_ts_file_line(self.Qs_ts_unit, self.time, self.Qs, IDs)

    #   save_pixel_values()
    #-------------------------------------------------------------------
#-----------------------------------------------------------------------
def Get_Timestep(cmin, nmin, dx, Rmax, Amax, Smax):

    #------------------------------------------------------------
    # Notes: The Courant condition: v < dx/dt is used together
    #        with the following equations to compute a stable
    #        timestep:

    #        (1)  Q = R * A = v * w * d
    #        (2)  v = d^(2/3) * S^(1/2) / n
    #        (2b) d = (n v )^(3/2) * S^(-3/4)

    #        Combining (1) and (2) we get:
    #
    #        (3) v  = (R * A / w)^(2/5) * S^(3/10) * n^(-3/5)
    #        (4) w  = c * dx     (c <= 1)
    #        (5) v  < dx / dt

    #        Combining these and solving for dt we get:

    #        (6) dt < [c^(2/5) * n^(3/5) * dx^(7/5)] /
    #                 [(R * A)^(2/5) * S^(3/10)]

    #        Use cmin, nmin, dx_min, Rmax, Amax, Smax.
    #------------------------------------------------------------
    numer = (cmin ** float64(0.4)) * (nmin ** float64(0.6)) * dx ** float64(1.4)
    denom = (Rmax * Amax) ** float64(0.4) * Smax ** float64(0.3)
    dt = (numer / denom)
    
    return dt
    
#   Get_Timestep()
#-----------------------------------------------------------------------
def Stable_Timestep(A, S, dx=None, dy=None, dt=None, R=None,
                    theta=None, m=None, n=None, k=None):

    #--------------------------------------------------------
    # Notes: This routine is based on a similarity idea for
    #        finding a stable timestep using the fact that
    #        the model was stable for a previous set of
    #        parameters.  Recall that:  Qs = K Q^m S^n,
    #        Q = R A^theta, qs = Q/dw

    #        K R^m A^(theta * m) S^n dt / (ds * dw) = const

    #        We also assume that dx=dy and:
    #            ds0/ds = dw0/dw = dx0/dx = dy0/dy
    #--------------------------------------------------------

    #-----------------------
    # Stable parameter set
    #-----------------------
    n_params = 2
    dx0 = float32(40.0)  #(meters)
    dy0 = float32(40.0)  #(meters)
    dt0 = float32(10.0)  #(years)
    R0 = float32(1.0)   #(m/year)
    m0 = float32(1.0)
    n0 = float32(1.0)
    k0 = float32(0.01)
    theta0 = float32(1.0)
    P0 = float32(100.0)  #(Not known as well, but max(A*S).)
    
    #-------------------
    # Keyword defaults
    #-------------------
    if (dx in [0,None]):    
        dx = dx0
    if (dy in [0,None]):    
        dy = dy0
    if (dt in [0,None]):    
        dt = dt0
    if (R in [0,None]):    
        R = R0
    if (m in [0,None]):    
        m = m0
    if (n in [0,None]):    
        n = n0
    if (k in [0,None]):    
        k = k0
    if (theta in [0,None]):    
        theta = theta0
    
    #------------------------------------
    # Get max value of the A-S function
    #------------------------------------
    grid = (A ** (theta * m)) * (S ** n)
    P = nanmax(grid)
    
    #---------------
    # New timestep
    #---------------
    dt = (dx / dx0) * (dx / dx0) * (k0 / k) * ((R0 ** m0) / (R ** m))
    dt = dt * (P0 / P) * dt0
    
    return int16(dt)
    
#   Stable_Timestep()
#-----------------------------------------------------------------------
