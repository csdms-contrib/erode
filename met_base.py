
## Copyright (c) 2001-2009, Scott D. Peckham
## January 2009  (converted from IDL)
## August 2009

#########################################################
#  NB!  self.PRECIP_ONLY toggle could be used to avoid
#       updating other met vars when not necessary.

#       "met process" now has its own timestep.
#########################################################   

#-----------------------------------------------------------------------
#  NOTES:  This file defines a "base class" for meteorology
#          components as well as any functions used by most or
#          all meteorology methods.  The methods of this class
#          should be over-ridden as necessary for different
#          methods of modeling meteorology.
#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class met_component     (inherits from CSDMS_base.py)
#      get_input_items()     # (not ready yet)
#      get_output_items()
#      set_constants()
#      ---------------------
#      initialize()
#      update()
#      finalize()
#      read_config_file()
#      ---------------------
#      get_cca_ports()           # (in CSDMS_base.py)
#      embed_child_components()
#      add_child_ports()
#      initialize_ports()
#      release_cca_ports()       # (in CSDMS_base.py)
#      -----------------------------------
#      update_P_integral()
#      update_P_max()
#      update_snow_vars()
#      -----------------------------------
#      update_richardson_number()
#      update_bulk_exchange_coeff()
#      update_sensible_heat_flux()
#      update_saturation_vapor_pressure()
#      update_vapor_pressure()
#      update_latent_heat_flux()
#      update_conduction_heat_flux()
#      update_advection_heat_flux()

#      update_net_shortwave_radiation()     # (not written yet, read as input)
#      update_net_longwave_radiation()      # (not written yet, read as input)
#      update_net_energy_flux()             # ("Q_sum")
#      -----------------------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      -----------------------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
#-----------------------------------------------------------------------

from numpy import *
import numpy

import CSDMS_base
import model_input
import tf_utils

from save_load import *  # (used by read_config_file())
from model_output import *

#-----------------------------------------------------------------------
def unit_test():

    directory   = tf_utils.TF_Test_Directory()
    data_prefix = tf_utils.TF_Test_Data_Prefix()
    case_prefix = tf_utils.TF_Test_Case_Prefix()
  
    m = met_component()
    m.CCA = False

    #--------------------------------------
    # Call initialize() and update() once
    #--------------------------------------
##    m.initialize(directory=directory,
##                 data_prefix=data_prefix,
##                 case_prefix=case_prefix, mode="main")
##    time_sec = float64(0)
##    m.update(time_sec)

    m.run_model(directory=directory,
                data_prefix=data_prefix,
                case_prefix=case_prefix)  
 
##    print 'nx          =', m.nx
##    print 'ny          =', m.ny
##    print 'config_file =', m.config_file
    
#   unit_test()
#-----------------------------------------------------------------------
class met_component(CSDMS_base.CSDMS_component):

        #-------------------------------------------------------------
        # Notes: **************************************************
        #        Do we ever need to distinguish between a surface
        #        temperature and snow temperature (in the snow) ?
        #        Recall that a separate T_soil_x variable is used
        #        to compute Qc.
        #        **************************************************

        #        Vapor_Pressure function is defined in Qnet_file.py
        #        and accepts temperature argument, in Celsius, and
        #        relative humidity.  It returns vapor pressure in
        #        units of kPa.

        #        Cp_snow is from NCAR CSM Flux Coupler web page

        #        uz = wind speed, z = reference height above ground

        #        rho_H2O is currently not adjustable with GUI.
        #------------------------------------------------------------

    #-------------------------------------------------------------------
    def set_constants(self):

        #---------------------------------
        # Define some physical constants
        #---------------------------------
        self.g        = float64(9.81)    # [m/s^2, gravity]
        self.kappa    = float64(0.408)   # [none, von Karman]
        self.rho_H2O  = float64(1000)
        self.rho_air  = float64(1.2614)
        self.Cp_air   = float64(1005.7)
        self.Lv       = float64(2500000)   #[Joules/kg, Latent heat of vaporiz.)
        self.Lf       = float64(334000)    #[J/kg = W*s/kg, Latent heat of fusion)

        #----------------------------------------
        # Constants related to precip (9/24/09)
        #----------------------------------------
        self.PRECIP_ONLY = False  #######  For efficiency later ? #####
        self.mmph_to_mps = (float64(1) / float64(3600000))
        self.mps_to_mmph = float64(3600000)
        self.forever     = float64(999999999)  # [minutes]
        
        #------------------------------------------------
        # Only needed for method 1, where all rates and
        # durations are read as 1D arrays from GUI.
        # Method 1 may be removed in a future version.
        #------------------------------------------------
##        self.method1_rates     = None
##        self.method1_durations = None
##        self.method1_n_rates   = 0
        
    #   set_constants()             
    #-------------------------------------------------------------------
    def initialize(self, directory=None,
                   data_prefix=None,
                   case_prefix=None, mode="module"):

        self.status = 'initializing'
        self.mode   = mode
        print 'Meteorology component: Initializing...'
        self.set_directory(directory, data_prefix, case_prefix)

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()       # (12/6/09)
        self.read_config_file()
        
        #------------------------------------------------------
        # NB! "Sample steps" must be defined before we return
        #     Check all other process modules.
        #------------------------------------------------------
        if (self.method == 0):
            self.P      = float64(0)
            self.e_air  = float64(0)   ## NEW SECTION
            self.e_surf = float64(0)
            self.Qn_SW  = float64(0)
            self.Qn_LW  = float64(0)
            self.Q_sum  = float64(0)
            self.DONE   = True
            self.status = 'initialized'
            return

        ### self.P = float64(0)
        
        self.read_grid_info()
        
        #----------------------------------------
        # Initialize all time-related variables
        #----------------------------------------
        self.initialize_time_vars()

        #-----------------------------------------------
        # Read from files as needed to initialize vars 
        #-----------------------------------------------
        self.open_input_files()
        self.read_input_files()  # (initializes P)

        #------------------
        # Initialize vars
        #-------------------------------------------------
        # Initialize max precip rate with the first rate
        #------------------------------------------------
        # Note: Need this here because rate may be
        #       zero at the end of update_precip_rate()
        #------------------------------------------------
        self.P_max = self.P.max()    # (after read_input_files)
        self.vol_P = float64(0)      # (for mass balance)
        self.Q_sum = float64(0)      # (for energy_balance snow + ET)
        
        self.initialize_required_components(mode)
        
        #-------------------------------
        # Save these for writing output
        #-------------------------------
        self.store_outlet_IDs()
        
        self.open_output_files()
        
        self.status = 'initialized'  # (OpenMI 2.0 convention)

    #   initialize()
    #-------------------------------------------------------------------
    def update(self, time_seconds=None):

        #----------------------------------------------------------
        # Note: The read_input_files() method is first called by
        #       the initialize() method.  Then, the update()
        #       method is called one or more times, and it calls
        #       other update_*() methods to compute additional
        #       variables using input data that was last read.
        #       Based on this pattern, read_input_files() should
        #       be called at end of update() method as done here.
        #       If the input files don't contain any additional
        #       data, the last data read persists by default.
        #----------------------------------------------------------
        if (self.method == 0):
            return
        self.status = 'updating'  # (OpenMI 2.0 convention)
        
        ###########################################################
        # Do we really need the tests shown here before calling
        # update routines, or can they handle this themselves?
        ###########################################################
        
##        #---------------------------------------------
##        # Compute new e_air from new T_air and RH
##        # if var type is Time Series of Grid Stack
##        #---------------------------------------------
##        if (self.T_air_type in [1,3]) or \
##           (self.RH_type in [1,3]):
##            self.update_saturation_vapor_pressure(MBAR=True)
##            self.update_vapor_pressure(MBAR=True)
##
##        #---------------------------------------------
##        # Compute new e_surf [kPa] from new T_surf
##        #---------------------------------------------
##        if (self.T_surf_type in [1,3]):
##            #----------------------------------------------------------
##            # Uses Saturation_Vapor_Pressure() from energy_balance.py
##            #----------------------------------------------------------
##            self.e_surf = Saturation_Vapor_Pressure(self.T_surf)  #[kPa]
##            # self.update_saturation_vapor_pressure_at_surface(MBAR=False)

        #-------------------------------------------
        # Update computed values related to precip
        #-------------------------------------------
        self.update_P_integral()
        self.update_P_max()
        self.update_snow_vars()   # (snow depth due to snowfall)
        
        #-------------------------
        # Update computed values
        #-------------------------
        if not(self.PRECIP_ONLY):
            self.update_richardson_number()
            self.update_bulk_exchange_coeff()
            self.update_sensible_heat_flux()
            self.update_saturation_vapor_pressure(MBAR=True)
            self.update_saturation_vapor_pressure(MBAR=True, SURFACE=True)  ##########
            self.update_vapor_pressure(MBAR=True)
            self.update_vapor_pressure(MBAR=True, SURFACE=True)   ########
            self.update_latent_heat_flux()          # (uses e_air and e_surf)
            self.update_conduction_heat_flux()
            self.update_advection_heat_flux()
            self.update_net_energy_flux()
            # self.update_net_shortwave_radiation()
            # self.update_net_longwave_radiation()

        #----------------------------------------
        # Read next met vars from input files ?
        #-------------------------------------------
        # Note that read_input_files() is called
        # by initialize() and these values must be
        # used for "update" calls before reading
        # new ones.
        #-------------------------------------------
        if (self.time_index > 0):
            self.read_input_files()
      
        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        self.write_output_files(time_seconds)

        #------------------------
        # Update internal clock
        #------------------------
        self.update_time()
        self.status = 'updated'  # (OpenMI)
        
    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI)
        if (self.method != 0):
            self.close_input_files()   ##  TopoFlow input "data streams"
            self.close_output_files()
        self.status = 'finalized'  # (OpenMI)
        print 'Meteorology component: Finished.'
        
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
        print 'Meteorology component: Reading config file...'
        self.config_file = (self.directory +
                            self.case_prefix + '_meteorology.cfg')
        file_unit = open(self.config_file, 'r')

        #-----------------------------
        # Skip over the header lines
        #-----------------------------
        for k in xrange(1, 5):
            line = file_unit.readline()
            
        #-------------------------------
        # Read the meteorological vars
        #-------------------------------
        method      = Read_Vars(file_unit, data_type='BYTE')
        method_name = Read_Vars(file_unit, data_type='STRING')
        dt_type, dt = Read_Vars(file_unit, True)
        #---------------------------------------------------
        self.method      = method
        self.method_name = method_name
        self.dt          = float64(dt)
        #---------------------------------------------------
        rho_H2O_type, rho_H2O = Read_Vars(file_unit, True)
        rho_air_type, rho_air = Read_Vars(file_unit, True)
        Cp_air_type, Cp_air   = Read_Vars(file_unit, True)
        #---------------------------------------------------
        self.rho_H2O = float64(rho_H2O)
        self.rho_air = float64(rho_air)
        self.Cp_air  = float64(Cp_air)
        #--------------------------------------------------
        P_type,     P       = Read_Vars(file_unit, True) ############
        T_air_type, T_air   = Read_Vars(file_unit, True)
        T_surf_type, T_surf = Read_Vars(file_unit, True)
        RH_type, RH         = Read_Vars(file_unit, True)
        p0_type, p0         = Read_Vars(file_unit, True)
        uz_type, uz         = Read_Vars(file_unit, True)
        z_type, z           = Read_Vars(file_unit, True)
        z0_air_type, z0_air = Read_Vars(file_unit, True)
        Qn_SW_type, Qn_SW   = Read_Vars(file_unit, True)
        Qn_LW_type, Qn_LW   = Read_Vars(file_unit, True)
        #--------------------------------------------------------
        P, P_file         = Load_Var(P, P_type)
        self.P            = P
        self.P_type       = Type_Code(P_type)
        self.P_file       = P_file
        
        T_air, T_air_file = Load_Var(T_air, T_air_type)
        self.T_air        = T_air
        self.T_air_type   = Type_Code(T_air_type)
        self.T_air_file   = T_air_file
        
        T_surf, T_surf_file = Load_Var(T_surf, T_surf_type)
        self.T_surf       = T_surf
        self.T_surf_type  = Type_Code(T_surf_type)
        self.T_surf_file  = T_surf_file
        
        RH, RH_file       = Load_Var(RH, RH_type)
        self.RH           = RH
        self.RH_type      = Type_Code(RH_type)
        self.RH_file      = RH_file
        
        p0, p0_file       = Load_Var(p0, p0_type)
        self.p0           = p0
        self.p0_type      = Type_Code(p0_type)
        self.p0_file      = p0_file
        
        uz, uz_file       = Load_Var(uz, uz_type)
        self.uz           = uz
        self.uz_type      = Type_Code(uz_type)
        self.uz_file      = uz_file
        
        z, z_file         = Load_Var(z, z_type)
        self.z            = z
        self.z_type       = Type_Code(z_type)
        self.z_file       = z_file
        
        z0_air, z0_air_file = Load_Var(z0_air, z0_air_type)
        self.z0_air       = z0_air
        self.z0_air_type  = Type_Code(z0_air_type)
        self.z0_air_file  = z0_air_file
        #--------------------------------------------------------
        Qn_SW, Qn_SW_file = Load_Var(Qn_SW, Qn_SW_type)
        self.Qn_SW        = Qn_SW
        self.Qn_SW_type   = Type_Code(Qn_SW_type)
        self.Qn_SW_file   = Qn_SW_file
        Qn_LW, Qn_LW_file = Load_Var(Qn_LW, Qn_LW_type)
        self.Qn_LW        = Qn_LW
        self.Qn_LW_type   = Type_Code(Qn_LW_type)
        self.Qn_LW_file   = Qn_LW_file
        
        #--------------------------------------------------------   
        # Everything below this line is new on August 5, 2009.
        # Completely removed from the snow process component.
        #--------------------------------------------------------
        # Later, add option to save a computed Qnet_SW and _LW.
        #-------------------------------------------------------------------------
        dum_str, save_grid_dt = Read_Vars(file_unit, True)   ## NEW
        save_ea_grids, ea_gs_file = Read_Vars(file_unit, True, data_type='BYTE')
        save_es_grids, es_gs_file = Read_Vars(file_unit, True, data_type='BYTE')
        #-------------------------------------------------------------------------
        self.save_grid_dt  = float64(save_grid_dt)
        self.SAVE_EA_GRIDS = save_ea_grids
        self.SAVE_ES_GRIDS = save_es_grids
        #------------------------------------
        self.ea_gs_file = ea_gs_file
        self.es_gs_file = es_gs_file
        #---------------------------------------------------------------------------        
        dum_str, save_pixels_dt     = Read_Vars(file_unit, True)
        save_ea_pixels, ea_ts_file = Read_Vars(file_unit, True, data_type='BYTE')
        save_es_pixels, es_ts_file = Read_Vars(file_unit, True, data_type='BYTE')
        #---------------------------------------------------------------------------
        self.save_pixels_dt = float64(save_pixels_dt)       
        self.SAVE_EA_PIXELS = save_ea_pixels
        self.SAVE_ES_PIXELS = save_es_pixels
        #-------------------------------------------- 
        self.ea_ts_file = ea_ts_file
        self.es_ts_file = es_ts_file
        
        #-----------------------
        # Close the input file
        #-----------------------
        file_unit.close()

        #---------------------------------------------------------
        # Make sure that all "save_dts" are larger or equal to
        # the specified process dt.  There is no point in saving
        # results more often than they change.
        # Issue a message to this effect if any are smaller ??
        #---------------------------------------------------------
        self.save_grid_dt   = maximum(self.save_grid_dt,    self.dt)
        self.save_pixels_dt = maximum(self.save_pixels_dt,  self.dt)
        
    #   read_config_file()
    #-------------------------------------------------------------------
    def embed_child_components(self):

        #------------------------------------------------
        # Note: Don't call this if (self.CCA == True).
        #------------------------------------------------
        # Instantiate and embed "process components"
        # in the place of the CCA ports.
        #------------------------------------------------
        # But how do we choose a given process method
        # such as "kinematic wave" ??
        #------------------------------------------------
        import basins
        import snow_degree_day
        
        self.bp = basins.basins_component()
        self.sp = snow_degree_day.snow_component()
        ###  self.pp = precip_base.precip_component()
        
    #   embed_child_components()
    #-------------------------------------------------------------------
    def add_child_ports(self):

        self.add_child_port('sp', 'mp', SELF=True)
        self.add_child_port('sp', 'bp')
        
    #   add_child_ports()        
    #-------------------------------------------------------------------
    def initialize_ports(self):
        
        #-----------------------------------------------
        # Initialize the process objects/components
        # This is also where output files are opened.
        #-----------------------------------------------
        # bp must be initialized first, since other inits
        # use the outlet_ID, etc.
        #--------------------------------------------------        
        if (self.bp.get_status() != 'initialized'):        # basin vars
            self.bp.initialize(directory=self.directory,
                               data_prefix=self.data_prefix,
                               case_prefix=self.case_prefix)
            
        if (self.sp.get_status() != 'initialized'):        # snow vars         
            self.sp.initialize(directory=self.directory,
                               data_prefix=self.data_prefix,
                               case_prefix=self.case_prefix)

##        if (self.pp.get_status() != 'initialized'):        # precip vars         
##            self.pp.initialize(directory=self.directory,
##                               data_prefix=self.data_prefix,
##                               case_prefix=self.case_prefix)
            
    #   initialize_ports()
    #-------------------------------------------------------------------
    def update_P_integral(self):

        #---------------------------------------------------
        # Notes: This can be used for mass balance checks,
        #        such as now done by update_mass_totals()
        #        in topoflow.py.  The "dt" here should be
        #        TopoFlow's "main dt" vs. the process dt.
        
        #        dV[i] = P[i] * da[i] * dt, dV = sum(dV[i])
        #---------------------------------------------------
        
        #------------------------------------------------
        # Update mass total for P, sum over all pixels
        #------------------------------------------------   
        volume = double(self.P * self.da * self.dt)  # [m^3]
        if (size(volume) == 1):
            self.vol_P += (volume * self.rti.n_pixels)
        else:
            self.vol_P += sum(volume)
        
    #   update_P_integral()
    #-------------------------------------------------------------------
    def update_P_max(self):

        #-----------------------------------------
        # Save the maximum precip. rate in [m/s]
        #-----------------------------------------
        self.P_max = maximum(self.P_max, self.P.max())

    #   update_P_max()   
    #-------------------------------------------------------------------
    def update_snow_vars(self):

        #---------------------------------------------
        # If P or T_air is a grid, then we must have
        # h_swe and h_snow be grids.  This is set
        # up at start of Route_Flow.
        #---------------------------------------------
        h_swe    = self.get_port_data('h_swe',  self.sp)
        h_snow   = self.get_port_data('h_snow', self.sp)
        ## print 'MADE IT PAST get_port_data() calls.'

        #--------------------------------------------------        
        rho_H2O  = self.sp.get_scalar_double('rho_H2O')
        rho_snow = self.sp.get_scalar_double('rho_snow')
        ### main_dt  = self.cp.get_scalar_double('dt')
        ## print 'MADE IT PAST get_scalar_double() calls.'
        
        #----------------------------------------------
        # Update snow water equivalent and snow depth
        #----------------------------------------------
        # (3/14/07) New method that works regardless
        # of whether P and T are scalars or grids.
        #---------------------------------------------
        # (8/20/09) Precip process now has constant
        # timestep, dt, which should be used here.
        #---------------------------------------------        
        dh_swe  = (self.P * (self.T_air <= 0)) * self.dt
        h_swe  += dh_swe
        ratio   = (rho_H2O / rho_snow)
        h_snow  = h_swe * ratio
        ## print 'MADE IT PAST var updates.'
        
        #----------------------------------------------------
        # Use "set_values" call to update h_snow and h_swe?
        #----------------------------------------------------
##        print 'shape(h_snow) =', shape(h_snow)
##        print 'shape(h_swe)  =', shape(h_swe)
##        print 'shape(P)      =', shape(self.P)
##        print 'shape(T_air)  =', shape(self.T_air)
##        print 'T_air         =', self.T_air
        
        self.set_port_data('h_snow', h_snow, self.sp)
        self.set_port_data('h_swe',  h_swe,  self.sp)
        ## print 'MADE IT PAST set_port_data() calls.'
        
        #---------------------------------------------
        # Where precip falls as snow and gets saved
        # to SWE, P must be set to zero so it can't
        # go on to produce runoff.
        #-----------------------------------------------
        # A special situation occurs if P is a scalar
        # and T is a grid.  P must then be converted
        # to a grid in order to be able to reflect the
        # fact that water was stored as snow in some
        # pixels and not at others.  This happens
        # automatically as written here. (3/14/07)
        # If P and T are both scalars, P stays scalar.
        ################################################
        # NB!  Should still work in Python, because:
        #   (1.5 * True = 1.5, 1.5 * False = 0)
        ################################################
        # NB!  This could be slowing things down a lot
        #      because it happens every timestep, not
        #      just when precip rate is updated.  Is
        #      there some way to test out of this ??
        ################################################
        ################################################
        self.P = (self.P * (dh_swe == 0))

        #--------------
        # For testing
        #--------------
        # nT = size(mv.T_air)
        # nP = size(self.P)
        # if (nT == 1):
        #     tstr = str(mv.T_air)
        #     print '   T_air = ' + tstr
        # if (nP == 1):
        #     pstr = str(self.P)
        #     print '   P     = ' + pstr
        
        #-------------------------------------
        # Message about rainfall or snowfall
        #-------------------------------------
        # Note that rain and snow may fall
        # simultaneously at different pixels
        #-------------------------------------
        # Reporting this slows model too much
        # but can be used for testing
        #-------------------------------------
        # w1 = where(dh_swe > 0)
        # n1 = size(w1[0])
        # if (n1 > 0): print '   Snow is falling...'
        #--------------------------------------------------
        # dh_rain = (P * (mv.T_air > 0))
        # w2 = where(dh_rain > 0)
        # n2 = size(w2[0])
        # if (n2 > 0): print '   Rain is falling...'

    #   update_snow_vars() 
    #-------------------------------------------------------------------
    def update_richardson_number(self):

        #---------------------------------------------------------------
        # Notes: Other definitions are possible, such as the one given
        #        by Dingman (2002, p. 599).  However, this one is the
        #        one given by Zhang et al. (2000) and is meant for use
        #        with the stability criterion also given there.
        #---------------------------------------------------------------
        top     = self.g * self.z * (self.T_air - self.T_surf)
        bot     = (self.uz)**2.0 * (self.T_air + float64(273.15))
        self.Ri = (top / bot)

    #   update_richardson_number()
    #-------------------------------------------------------------------        
    def update_bulk_exchange_coeff(self):

        #----------------------------------------------------------------
        # Notes: Dn       = bulk exchange coeff for the conditions of
        #                   neutral atmospheric stability [m/s]
        #        Dh       = bulk exchange coeff for heat  [m/s]
        #        De       = bulk exchange coeff for vapor [m/s]
        #        h_snow   = snow depth [m]
        #        z0_air   = surface roughness length scale [m]
        #                   (includes vegetation not covered by snow)
        #        z        = height that has wind speed uz [m]
        #        uz       = wind speed at height z [m/s]
        #        kappa    = 0.408 = von Karman's constant [unitless]
        #        RI       = Richardson's number (see function)
        #----------------------------------------------------------------
        h_snow = self.get_port_data('h_snow', self.sp)
        
        #---------------------------------------------------
        # Compute bulk exchange coeffs (neutral stability)
        # using the logarithm "law of the wall".
        #---------------------------------------------------
        arg = self.kappa / log((self.z - h_snow) / self.z0_air)
        Dn  = self.uz * (arg)**2.0
        
        #-----------------------------------------------
        # NB! Dn could be a scalar or a grid, so this
        #     must be written to handle both cases.
        #     Note that WHERE can be used on a scalar:
        
        #     IDL> a = 1
        #     IDL> print, size(a)
        #     IDL> w = where(a ge 1, nw)
        #     IDL> print, nw
        #     IDL> a[w] = 2
        #     IDL> print, a
        #     IDL> print, size(a)
        #-----------------------------------------------

        ###########################################################
        #  NB!  If T_air and T_surf are both scalars, then next
        #       few lines won't work because we can't index the
        #       resulting empty "w" (even if T_air == T_surf).
        ###########################################################
##        w  = where(self.T_air != self.T_surf)
##        nw = size(w[0])
##        ## nw = size(w,0)  # (doesn't work if 2 equal scalars)
        #----------------------------------------------------------
        T_AIR_SCALAR  = (rank(self.T_air)  == 0)
        T_SURF_SCALAR = (rank(self.T_surf) == 0)
        if (T_AIR_SCALAR and T_SURF_SCALAR):
            if (self.T_air == self.T_surf):  nw=1
            else: nw=0      
        else:
            w  = where(self.T_air != self.T_surf)
            nw = size(w[0])
        #----------------------------------------------------------            
        if (nw == 0):
            self.Dn = Dn
            self.Dh = Dn
            self.De = Dn
            return
        
        #-------------------------------------
        # One or more pixels are not neutral
        # so make a correction using RI
        #---------------------------------------------
        # NB!  RI could be a grid when Dn is a
        # scalar, and this will change Dn to a grid.
        #---------------------------------------------
        # Ri = Richardson_Number(z, uz, T_air, T_surf)
        #--------------------------------------------
        # Before 12/21/07.  Has bug if RI is a grid
        #--------------------------------------------
        # w_stable = where(*T_air gt *T_surf, n_stable)
        # if (n_stable ne 0) then begin
        #     Dn[w_stable] = Dn[w_stable]/(1d + (10d * RI))
        # endif
        # w_unstable = where(*T_air lt *T_surf, n_unstable)
        # if (n_unstable ne 0) then begin
        #----------------------------------------------
        # Multiplication and substraction vs. opposites
        # for the stable case.  Zhang et al. (2000)
        # Hopefully not just a typo.
        #----------------------------------------------
        #    Dn[w_unstable] = Dn[w_unstable]*(1d - (10d * self.Ri))
        # endif
        
        #-----------------
        # After 12/21/07
        #------------------------------------------------------------
        # If T_air, T_surf or uz is a grid, then Ri will be a grid.
        # This version makes only one call to WHERE, so its faster.
        #------------------------------------------------------------
        # Multiplication and substraction vs. opposites for the
        # stable case (Zhang et al., 2000); hopefully not a typo.
        # It plots as a smooth curve through Ri=0.
        #------------------------------------------------------------
        nD = size(Dn)
        nR = size(self.Ri)
        if (nR > 1):    
            #--------------------------
            # Case where RI is a grid
            #--------------------------
            ws = where(self.Ri > 0)
            ns = size(ws[0])
            wu = where(invert(self.Ri > 0))
            nu = size(wu[0])
            if (nD == 1):    
                #******************************************
                # Convert Dn to a grid here or somewhere
                # Should stop with an error message
                #******************************************
                dum = int16(0)
            if (ns != 0):    
                Dn[ws] = Dn[ws] / (float64(1) + (float64(10) * self.Ri[ws]))
            if (nu != 0):    
                Dn[wu] = Dn[wu] * (float64(1) - (float64(10) * self.Ri[wu]))
        else:    
            #----------------------------
            # Case where Ri is a scalar
            #--------------------------------
            # Works if Dn is grid or scalar
            #--------------------------------
            if (self.Ri > 0):    
                Dn = Dn / (float64(1) + (float64(10) * self.Ri))
            else:    
                Dn = Dn * (float64(1) - (float64(10) * self.Ri))

        #----------------------------------------------------
        # NB! We currently assume that these are all equal.
        #----------------------------------------------------
        self.Dn = Dn
        self.Dh = Dn
        self.De = Dn
        
    #   update_bulk_exchange_coeff()
    #-------------------------------------------------------------------
    def update_sensible_heat_flux(self):

        #--------------------------------------------------------
        # Notes: All the Q's have units of W/m^2 = J/(m^2 s).
        #        Dh is returned by Bulk_Exchange_Coeff function
        #        and is not a pointer.
        #--------------------------------------------------------
      
        #---------------------
        # Physical constants
        #---------------------
        # rho_air = 1.225d   ;[kg/m^3, at sea-level]
        # Cp_air  = 1005.7   ;[J/kg/deg_C]
        
        #-----------------------------
        # Compute sensible heat flux
        #-----------------------------
        delta_T = (self.T_air - self.T_surf)
        self.Qh = (self.rho_air * self.Cp_air) * self.Dh * delta_T

    #   update_sensible_heat_flux()
    #-------------------------------------------------------------------
    def update_saturation_vapor_pressure(self, MBAR=False,
                                         SATTERLUND=False,
                                         SURFACE=False):

        #----------------------------------------------------------------
        #Notes:  Saturation vapor pressure is a function of temperature.
        #        T is temperature in Celsius.  By default, the method
        #        of Brutsaert (1975) is used.  However, the SATTERLUND
        #        keyword is set then the method of Satterlund (1979) is
        #        used.  When plotted, they look almost identical.  See
        #        the Compare_em_air_Method routine in Qnet_file.pro.
        #        Dingman (2002) uses the Brutsaert method.
        #        Liston (1995, EnBal) uses the Satterlund method.

        #        By default, the result is returned with units of kPa.
        #        Set the MBAR keyword for units of millibars.
        #        100 kPa = 1 bar = 1000 mbars
        #                => 1 kPa = 10 mbars
        #----------------------------------------------------------------
        #NB!     Here, 237.3 is correct, and not a misprint of 273.2.
        #        See footnote on p. 586 in Dingman (Appendix D).
        #----------------------------------------------------------------
        if (SURFACE):
            T = self.T_surf
        else:
            T = self.T_air
        
        if not(SATTERLUND):    
            #------------------------------
            # Use Brutsaert (1975) method
            #------------------------------
            term1 = (float64(17.3) * T) / (T + float32(237.3))
            e_sat = float64(0.611) * exp(term1)        #[kPa]
        else:    
            #-------------------------------
            # Use Satterlund (1979) method
            #-------------------------------
            e_sat = float64(10) ** (float64(11.4) - (float64(2353) / (T + float32(273.15))))   #[Pa]
            e_sat = (e_sat / float64(1000))   #[kPa]

        #-----------------------------------
        # Convert units from kPa to mbars?
        #-----------------------------------
        if (MBAR):    
            e_sat = (e_sat * float64(10))   #[mbar]

        if (SURFACE):
            self.e_sat_surf = e_sat
        else:
            self.e_sat_air  = e_sat

    #   update_saturation_vapor_pressure()
    #------------------------------------------------------------------- 
    def update_vapor_pressure(self, MBAR=False,
                              SATTERLUND=False,
                              SURFACE=False):

        #---------------------------------------------------
        # Notes: T is temperature in Celsius
        #        RH = relative humidity, in [0,1]
        #             by definition, it equals (e / e_sat)
        #        e has units of kPa.
        #---------------------------------------------------
        if (SURFACE):
            e_sat = self.e_sat_surf
        else:
            e_sat = self.e_sat_air
            
        e = (self.RH * e_sat)
        
        #-----------------------------------
        # Convert units from kPa to mbars?
        #-----------------------------------
        if (MBAR):    
            e = (e * float64(10))   #[mbar]

        if (SURFACE):
            self.e_surf = e
        else:
            self.e_air  = e
 
    #   update_vapor_pressure()
    #-------------------------------------------------------------------
    def update_latent_heat_flux(self):

        #---------------------------
        # Compute latent heat flux
        #----------------------------------------------------------
        # Notes: Pressure units cancel out.
        # According to Dingman (2002, p. 273), constant should be
        # 0.622 instead of 0.662 (Zhang et al., 2000).
        #----------------------------------------------------------
        # const = float64(0.622)
        const   = float64(0.662)
        factor  = (self.rho_air * self.Lv * self.De)
        delta_e = (self.e_air - self.e_surf)
        self.Qe = factor * (const / self.p0) * delta_e

    #   update_latent_heat_flux()
    #-------------------------------------------------------------------
    def update_conduction_heat_flux(self):

        #-----------------------------------------------------------------
        # Notes: The conduction heat flux from snow to soil for computing
        #        snowmelt energy, Qm, is close to zero.

        #        However, the conduction heat flux from surface and sub-
        #        surface for computing Qet is given by Fourier's Law,
        #        namely Qc = Ks(Tx - Ts)/x.

        #        All the Q's have units of W/m^2 = J/(m^2 s).
        #-----------------------------------------------------------------
        self.Qc = float64(0)

    #   update_conduction_heat_flux()
    #-------------------------------------------------------------------
    def update_advection_heat_flux(self):

        #------------------------------------------------------
        # Notes: All the Q's have units of W/m^2 = J/(m^2 s).
        #------------------------------------------------------
        self.Qa = float64(0)
        
    #   update_advection_heat_flux()
    #-------------------------------------------------------------------
    def update_net_shortwave_radiation(self):

        pass  # (not ready yet)
    
    #   update_net_shortwave_radiation()
    #-------------------------------------------------------------------
    def update_net_longwave_radiation(self):

        pass  # (not ready yet)

    #   update_net_longwave_radiation()
    #-------------------------------------------------------------------
    def update_net_energy_flux(self):

        #------------------------------------------------------
        # Notes: Q_sum is used by "snow_energy_balance.py".
        #------------------------------------------------------
        #        Qm    = energy used to melt snowpack (if > 0)
        #        Qn_SW = net shortwave radiation flux (solar)
        #        Qn_LW = net longwave radiation flux (air, surface)
        #        Qh    = sensible heat flux from turbulent convection
        #                between snow surface and air
        #        Qe    = latent heat flux from evaporation, sublimation,
        #                and condensation
        #        Qa    = energy advected by moving water (i.e. rainfall)
        #                (ARHYTHM assumes this to be negligible; Qa=0.)
        #        Qc    = energy flux via conduction from snow to soil
        #                (ARHYTHM assumes this to be negligible; Qc=0.)
        #        Ecc   = cold content of snowpack = amount of energy
        #                needed before snow can begin to melt [J/m^2]

        #        All Q's here have units of [W/m^2].
        #        Are they all treated as positive quantities ?

        #        rho_air  = density of air [kg/m^3]
        #        rho_snow = density of snow [kg/m^3]
        #        Cp_air   = specific heat of air [J/kg/deg_C]
        #        Cp_snow  = heat capacity of snow [J/kg/deg_C]
        #                 = ???????? = specific heat of snow
        #        Kh       = eddy diffusivity for heat [m^2/s]
        #        Ke       = eddy diffusivity for water vapor [m^2/s]
        #        Lv       = latent heat of vaporization [J/kg]
        #        Lf       = latent heat of fusion [J/kg]
        #        ------------------------------------------------------
        #        Dn       = bulk exchange coeff for the conditions of
        #                   neutral atmospheric stability [m/s]
        #        Dh       = bulk exchange coeff for heat
        #        De       = bulk exchange coeff for vapor
        #        ------------------------------------------------------
        #        T_air    = air temperature [deg_C]
        #        T_surf   = surface temperature [deg_C]
        #        T_snow   = average snow temperature [deg_C]
        #        RH       = relative humidity [none] (in [0,1])
        #        e_air    = air vapor pressure at height z [mbar]
        #        e_surf   = surface vapor pressure [mbar]
        #        ------------------------------------------------------
        #        h_snow   = snow depth [m]
        #        z        = height where wind speed is uz [m]
        #        uz       = wind speed at height z [m/s]
        #        P0       = atmospheric pressure [mbar]
        #        T0       = snow temperature when isothermal [deg_C]
        #                   (This is usually 0.)
        #        z0_air   = surface roughness length scale [m]
        #                   (includes vegetation not covered by snow)
        #                   (Values from page 1033: 0.0013, 0.02 [m])
        #        kappa    = von Karman's constant [unitless] = 0.41
        #        dt       = snowmelt timestep [seconds]
        #----------------------------------------------------------------
        self.Q_sum = self.Qn_SW + self.Qn_LW + self.Qh + \
                     self.Qe + self.Qa + self.Qc    #[W/m^2]

    #   update_net_energy_flux()   
    #-------------------------------------------------------------------  
    def open_input_files(self):

        self.P_unit      = model_input.open_file(self.P_type,      self.P_file)
        self.T_air_unit  = model_input.open_file(self.T_air_type,  self.T_air_file)
        self.T_surf_unit = model_input.open_file(self.T_surf_type, self.T_surf_file)
        self.RH_unit     = model_input.open_file(self.RH_type,     self.RH_file)
        self.p0_unit     = model_input.open_file(self.p0_type,     self.p0_file)
        self.uz_unit     = model_input.open_file(self.uz_type,     self.uz_file)
        self.z_unit      = model_input.open_file(self.z_type,      self.z_file)
        self.z0_air_unit = model_input.open_file(self.z0_air_type, self.z0_air_file)
        #--------------------------------------------------------------------------
        self.Qn_SW_unit  = model_input.open_file(self.Qn_SW_type,  self.Qn_SW_file)
        self.Qn_LW_unit  = model_input.open_file(self.Qn_LW_type,  self.Qn_LW_file)
        
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        rti = self.rti

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        P = model_input.read_next(self.P_unit, self.P_type, rti,
                               factor=self.mmph_to_mps)
        if (P != None):
            self.P = float32(P)
            if (self.time_index == 0):
                P_max  = P.max()
                print 'P_max =', P_max, ' (in read_input_files())'

            #--------------
            # For testing
            #--------------
            print 'min(P) =', P.min() * self.mps_to_mmph, ' [mmph]'
            print 'max(P) =', P.max() * self.mps_to_mmph, ' [mmph]'
        else:
            #------------------------------------
            # Out of data, set rate to zero
            # Precip is unique in this respect.
            #------------------------------------
            self.P = float32(0)
            ## print 'EOF, so P set to 0 by read_input_files().'  
        
        T_air = model_input.read_next(self.T_air_unit, self.T_air_type, rti)
        if (T_air != None): self.T_air = T_air

        T_surf = model_input.read_next(self.T_surf_unit, self.T_surf_type, rti)
        if (T_surf != None): self.T_surf = T_surf

        RH = model_input.read_next(self.RH_unit, self.RH_type, rti)
        if (RH != None): self.RH = RH

        p0 = model_input.read_next(self.p0_unit, self.p0_type, rti)
        if (p0 != None): self.p0 = p0

        uz = model_input.read_next(self.uz_unit, self.uz_type, rti)
        if (uz != None): self.uz = uz

        z = model_input.read_next(self.z_unit, self.z_type, rti)
        if (z != None): self.z = z

        z0_air = model_input.read_next(self.z0_air_unit, self.z0_air_type, rti)
        if (z0_air != None): self.z0_air = z0_air

        #-------------------------------------------------------------
        # These are currently treated as input data, but are usually
        # generated by functions in Qnet_file.py.  Later on, we'll
        # provide the option to compute them "on the fly" with new
        # functions called "update_net_shortwave_radiation()" and
        # "update_net_longwave_radiation()", called from update().
        #-------------------------------------------------------------        
        Qn_SW = model_input.read_next(self.Qn_SW_unit, self.Qn_SW_type, rti)
        if (Qn_SW != None): self.Qn_SW = Qn_SW

        Qn_LW = model_input.read_next(self.Qn_LW_unit, self.Qn_LW_type, rti)
        if (Qn_LW != None): self.Qn_LW = Qn_LW
         
    #   read_input_files()
    #-------------------------------------------------------------------  
    def close_input_files(self):

        if (self.P_file      != ''): self.P_unit.close()
        if (self.T_air_file  != ''): self.T_air_unit.close()
        if (self.T_surf_file != ''): self.T_surf_unit.close()
        if (self.RH_file     != ''): self.RH_unit.close()
        if (self.p0_file     != ''): self.p0_unit.close()
        if (self.uz_file     != ''): self.uz_unit.close()
        if (self.z_file      != ''): self.z_unit.close()
        if (self.z0_air_file != ''): self.z0_air_unit.close()
        #--------------------------------------------------------
        if (self.Qn_SW_file  != ''): self.Qn_SW_unit.close()        
        if (self.Qn_LW_file  != ''): self.Qn_LW_unit.close()
        
    #   close_input_files()
    #-------------------------------------------------------------------  
    def update_outfile_names(self):

        #--------------------------------------------------------
        # Notes:  Whenever the case_prefix changes (e.g. user
        #         changes it and clicks the Next button in the
        #         "run var" panel), all output filenames should
        #         be reset to new defaults.
        #--------------------------------------------------------
        self.ea_gs_file = (self.case_prefix + '_2D-ea.rts')
        self.es_gs_file = (self.case_prefix + '_2D-es.rts')
        #-----------------------------------------------------
        self.ea_ts_file = (self.case_prefix + '_0D-ea.txt')
        self.es_ts_file = (self.case_prefix + '_0D-es.txt')

    #   update_outfile_names()   
    #-------------------------------------------------------------------  
    def open_output_files(self):
        
        #----------------------------------
        # Open files to write grid stacks
        #----------------------------------
        if (self.SAVE_EA_GRIDS):    
            self.ea_gs_unit = open_new_gs_file(self.ea_gs_file, self.rti,
                                               grid_name='e_air',
                                               long_name='vapor_pressure_in_air',
                                               units_name='mbar')
            
        if (self.SAVE_ES_GRIDS):    
            self.es_gs_unit = open_new_gs_file(self.es_gs_file, self.rti,
                                               grid_name='e_surf',
                                               long_name='vapor_pressure_at_surface',
                                               units_name='mbar')

        #---------------------------------------
        # Open text files to write time series
        #---------------------------------------
        IDs = self.outlet_IDs 
        if (self.SAVE_EA_PIXELS):
            self.ea_ts_unit = open_new_ts_file( self.ea_ts_file, 'e_air', IDs)

        if (self.SAVE_ES_PIXELS):
            self.es_ts_unit = open_new_ts_file( self.es_ts_file, 'e_surf', IDs)

    #   open_output_files()
    #-------------------------------------------------------------------
    def write_output_files(self, time_seconds=None):

        #-----------------------------------------
        # Allows time to be passed from a caller
        #-----------------------------------------
        if (time_seconds is None):
            time_seconds = self.time_sec
        model_time = int(time_seconds)
        
        #----------------------------------------
        # Save computed values at sampled times
        #----------------------------------------
        if (model_time % int(self.save_grid_dt) == 0):
            self.save_grids()
        if (model_time % int(self.save_pixels_dt) == 0):
            self.save_pixel_values()

    #  write_output_files()
    #-------------------------------------------------------------------
    def close_output_files(self):
    
        if (self.SAVE_EA_GRIDS): self.ea_gs_unit.close()
        if (self.SAVE_ES_GRIDS): self.es_gs_unit.close()
        #---------------------------------------------------------
        if (self.SAVE_EA_PIXELS): self.ea_ts_unit.close()
        if (self.SAVE_ES_PIXELS): self.es_ts_unit.close()

    #   close_output_files()        
    #-------------------------------------------------------------------  
    def save_grids(self):
        
        nx = self.nx
        ny = self.ny
       
        if (self.SAVE_EA_GRIDS):
            save_as_grid_to_file(self.ea_gs_unit, self.e_air,
                                 'e_air', nx, ny)
            
        if (self.SAVE_ES_GRIDS):
            save_as_grid_to_file(self.es_gs_unit, self.e_surf,
                                 'e_surf', nx, ny)

    #   save_grids()            
    #-------------------------------------------------------------------  
    def save_pixel_values(self):

        IDs      = self.outlet_IDs
        time_min = self.time_min
        
        if (self.SAVE_EA_PIXELS):
            write_ts_file_line(self.ea_ts_unit, time_min, self.e_air, IDs)
            
        if (self.SAVE_ES_PIXELS):
            write_ts_file_line(self.es_ts_unit, time_min, self.e_surf, IDs)

    #   save_pixel_values()
    #-------------------------------------------------------------------


    
        
