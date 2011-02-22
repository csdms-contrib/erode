
## Copyright (c) 2001-2009, Scott D. Peckham
## August 2009

#-----------------------------------------------------------------------
#  Notes:  This file defines a "base class" for CSDMS "process"
#          components.  Many of these functions are implementations
#          of methods defined in "topoflow3.IRFPort.sidl".  They are 
#          therefore required to use TF components in a CCA framework.
#-----------------------------------------------------------------------
#
#  unit_test()
#
#  class CSDMS_component
#      __init__()
#      get_status()
#      is_scalar()
#      is_vector()
#      is_grid()
#      has_variable()
#      -----------------------------
#       Next 3 currently identical
#      -----------------------------
#      get_scalar_double()
#      get_vector_double()          ## (2/16/10
#      get_grid_double()
#      get_values_in_grid_double()  ## (2/17/10)
#      -----------------------------
#       Next 3 currently identical
#      -----------------------------
#      set_scalar_double()
#      set_vector_double()          ## (2/16/10)
#      set_grid_double()
#      set_values_in_grid_double()  ## (2/17/10)
#      ---------------------
#      get_scalar_long()
#      get_vector_long()
#      get_grid_long()
#      get_values_in_grid_long()    ## (2/17/10)
#      ---------------------
#      set_scalar_long()            ## (2/16/10)
#      set_vector_long()            ## (2/16/10)
#      set_grid_long()              ## (2/16/10)
#      set_values_in_grid_long()    ## (2/17/10)
#      ---------------------
#      get_input_items()
#      get_output_items()
#      ---------------------
#      read_gui_info_file()     ### 11/13/09
#      get_user_input()         ### 9/24/09, 11/13/09
#      load_user_input()        ### 10/1/09
#      ---------------------
#      go()
#      run_model()
#      ### read_config_file()   # (template)
#      ### initialize   # (template)
#      ### update()     # (template)
#      finalize()

#      -----------------------------------------------
#      These methods are not part of "IRF" interface
#      but are used by the TopoFlow components.
#      -----------------------------------------------
#      initialize_required_components()
#      get_cca_port()
#      get_cca_ports()
#      release_cca_ports()
#      ------------------------
#      add_child_port()
#      get_port_data()     # (rename to get_port_double ??)
#      set_port_data()     # (rename to set_port_double ??)
#      ------------------------
#      get_rank()
#      get_size()
#      ------------------------
#      set_directory()
#      read_grid_info()
#      store_outlet_IDs()
#      initialize_time_vars()
#      update_time()
#      print_time_and_value()
#      print_run_time()

#-----------------------------------------------------------------------

from numpy import *
import numpy

import sys, os, time

import cfg_files as cfg
import pixels
import rti_files
import tf_utils

#-----------------------------------------------------------------------
def unit_test():

    directory   = tf_utils.TF_Test_Directory()
    data_prefix = tf_utils.TF_Test_Data_Prefix()
    case_prefix = tf_utils.TF_Test_Case_Prefix()

    c = CSDMS_component()
    c.CCA = False
    print 'Instantiated component.'

    #--------------------------------------------
    # Check the "read_gui_info_file()" function
    #--------------------------------------------
##    gui_info_file = '/Applications/Erode/gui_info/Erode_GUI_Info.txt'
##    gui_info_file = '/data/progs/erode/3.0/gui_info/Erode_GUI_Info.txt'
    gui_info_file = '/Users/peckhams/Desktop/GC2D_GUI_Info.txt'
    var_names, labels, values, min_vals, max_vals, \
               desc, group_names, group_sizes = \
                   c.read_gui_info_file( gui_info_file )
    print 'var_names ='
    print var_names
    print 'labels ='
    print labels
    print 'values ='
    print values
    print 'min_vals ='
    print min_vals
    print 'max_vals ='
    print max_vals
    print 'desc ='
    print desc
    print 'group_names ='
    print group_names
    print 'group_sizes ='
    print group_sizes    
    return

    #-------------------------------------
    # Test the "print_run_time()" method
    #-------------------------------------
##    print ' '
##    print 'Testing "print_run_time()"...'
##    c.print_run_time(seconds=1.618)
##    c.print_run_time(seconds=60)
##    c.print_run_time(seconds=3600)
##    c.print_run_time(seconds=3600 * 24)
##    c.print_run_time(seconds=3600 * 24 * 365)
    
    #---------------------------
    # Test some of the methods
    #---------------------------
    c.a = numpy.float64(5)
    print 'c.a = numpy.float64(5)'
    print "c.is_scalar('a') =", c.is_scalar('a')
    print "c.is_grid('a')   =", c.is_grid('a')
    v1 = c.get_port_data('a', c)
    print "c.get_port_data('a',c) =", v1
    print ' '
    #-------------------------------------------------
    c.b = numpy.zeros((3,3), dtype='Float64')
    print "c.b = numpy.zeros((3,3), dtype='Float64')"
    print "c.is_scalar('b') =", c.is_scalar('b')
    print "c.is_grid('b')   =", c.is_grid('b')
    v2 = c.get_port_data('b', c)
    print "c.get_port_data('b',c) =", v2
    print ' '
    #-------------------------------------------------
    print "c.is_scalar('b[1]') =", c.is_scalar('b[1]')
    print "c.is_grid('b[1]')   =", c.is_grid('b[1]')
    v3 = c.get_port_data('b[1]', c)
    print "c.get_port_data('b[1]',c) =", v3
    print ' '
    #-------------------------------------------------    
    
    # This component has no initialize() method
##    c.initialize(directory=directory,
##                 data_prefix=data_prefix,
##                 case_prefix=case_prefix)
##    print 'nx =', c.nx
##    print 'ny =', c.ny

#   unit_test()
#-----------------------------------------------------------------------
class CSDMS_component:

    def __init__(self):

        self.CCA    = tf_utils.TF_Use_CCA()
        self.DEBUG  = True
        self.SILENT = True
        self.REPORT = False
        self.status = 'created'   # (OpenMI 2.0 conventions)
        self.USER_SET_VALUES = False
        
    #   __init__()
    #-------------------------------------------------------------------
    def get_status(self):

        #-----------------------------------------------------
        # Notes: Return component status as a string.  The
        #        possible return values are from OpenMI 2.0:
        #
        #           created, initializing, initialized,
        #           updating, updated, finalizing, finalized,
        #           failed (could add "stopped").
        #-----------------------------------------------------
        return self.status

    #   get_status()
    #-------------------------------------------------------------------
    def is_scalar(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #-------------------------------------------------      
        exec("n = numpy.rank(self." + var_name + ")")       
        return (n == 0)
    
    #   is_scalar()
    #-------------------------------------------------------------------
    def is_vector(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #------------------------------------------------     
        exec("n = numpy.rank(self." + var_name + ")")       
        return (n == 1)
    
    #   is_vector()
    #-------------------------------------------------------------------
    def is_grid(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #------------------------------------------------ 

        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------         
##        exec("type_str = str(type(self." + var_name + "))")
##        p1 = type_str.find("ndarray")
##        p2 = type_str.find("float")
##        if (p1 == -1) and (p2 == -1):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------        
##        if ("ndarray" not in type_str) and \
##           ("float" not in type_str):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #-------------------------------------------------------        
        exec("n = numpy.rank(self." + var_name + ")")
        return (n == 2)

    #   is_grid()
    #-------------------------------------------------------------------
    def has_variable(self, var_name):

        #------------------------------------------------------
        # If var_name includes square brackets for subscripts
        # remove them to get the actual variable name.
        #------------------------------------------------------
        bracket_pos = var_name.find('[')
        if (bracket_pos != -1):
            key = var_name[0:bracket_pos]
        else:
            key = var_name

        #---------------------------------------------------
        # Does current component have requested variable ?
        #---------------------------------------------------
        SILENT = True
        VARIABLE_FOUND = self.__dict__.has_key(key)
        if not(VARIABLE_FOUND) and not(SILENT):
            print 'WARNING: Component does not have the'
            print '         requested variable: ' + var_name
            print ' '
            
        return VARIABLE_FOUND
        
    #   has_variable()
    #-------------------------------------------------------------------        
    def get_scalar_double(self, var_name):

        #------------------------------------
        # Note: The next line doesn't work.
        #------------------------------------
        ## exec("return self." + var_name)

        #---------------------------------------------------
        # Does current component have requested variable ?
        #---------------------------------------------------
        # This is not used yet because possible impact on
        # performance has not be tested yet. (2/17/10)
        # If it does get used later, it will be added to
        # all of the "getters".
        #---------------------------------------------------        
##        if not(self.has_variable(var_name)):
##            return float64(0)
 
        try:
            exec("result = self." + var_name)
            return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_scalar_double().'
            print '    Returning 0.'
            return numpy.float64(0)

            ############## flush output here ??
        
    #   get_scalar_double()
    #-------------------------------------------------------------------
    def get_vector_double(self, var_name):

        #---------------------------------------------------------
        # Note: This function was causing a "segmentation fault
        #       in gui-backend.sh" error message when trying to
        #       run TopoFlow through the CMT (in CCA framework).
        #       Solution was to use numpy.array, as shown.
        #       (2/17/10)
        #---------------------------------------------------------
        try:
            exec("result = self." + var_name)
            return numpy.array(result, dtype='float64')
            #-------------------------
            # NB! This doesn't work.
            #-------------------------
            # return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_vector_double().'
            print '    Returning zeros.'
            return numpy.zeros([1], dtype='float64')
        
    #   get_vector_double()
    #-------------------------------------------------------------------
    def get_grid_double(self, var_name):

        try:
            exec("result = self." + var_name)
            return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_grid_double().'
            print '    Returning zeros.'
            return numpy.zeros([1,1], dtype='float64')
        
    #   get_grid_double()
    #-------------------------------------------------------------------
    def get_values_in_grid_double(self, var_name, IDs):

        #---------------------------------------------------------
        # Note: This function was causing a "segmentation fault
        #       in gui-backend.sh" error message when trying to
        #       run TopoFlow through the CMT (in CCA framework).
        #       Solution was to use numpy.array, as shown.
        #       (2/18/10)
        #---------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #---------------------------------------------------------
        try:
            exec("result = self." + var_name + '.flat[IDs]')
            return numpy.array(result, dtype='float64')
            ## return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_values_in_grid_double().'
            print '    Returning zeros.'
            return numpy.zeros(len(IDs), dtype='float64')
        
    #   get_values_in_grid_double()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def set_scalar_double(self, var_name, scalar):

        exec("self." + var_name + " = numpy.float64(scalar)")
    
    #   set_scalar_double()
    #-------------------------------------------------------------------
    def set_vector_double(self, var_name, vector):

        #--------------------------------------------------
        # Notes: First method here should be more robust.
        #        See Notes for get_vector_double().
        #--------------------------------------------------
        exec("self." + var_name + " = numpy.array(vector, dtype='float64')")

        #-----------------------------------
        # Original method (before 2/17/10)
        #-----------------------------------        
        # exec("self." + var_name + " = numpy.float64(vector)")
    
    #   set_vector_double()
    #-------------------------------------------------------------------
    def set_grid_double(self, var_name, grid):
         
        exec("self." + var_name + " = numpy.float64(grid)")
    
    #   set_grid_double()
    #-------------------------------------------------------------------
    def set_values_in_grid_double(self, var_name, IDs, values):

        #--------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #--------------------------------------------------------

        exec("self." + var_name + ".flat[IDs] = values")
        
        # exec("self." + var_name + ".flat[IDs] = numpy.float64(values)")
        
    #   set_values_in_grid_double()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def get_scalar_long(self, var_name):

        exec("result = numpy.int32(self." + var_name + ")")
        return result
    
    #   get_scalar_long()
    #-------------------------------------------------------------------
    def get_vector_long(self, var_name):

        #--------------------------------------------
        # Notes: See Notes for get_vector_double().
        #--------------------------------------------
        try:
            exec("result = self." + var_name)
            return numpy.array(result, dtype='int32')
            #-------------------------
            # NB! This doesn't work.
            #-------------------------
            # return numpy.int32(result)
        except:
            print 'ERROR in CSDMS_base.get_vector_long().'
            print '    Returning zeros.'
            return numpy.zeros([1], dtype='int32')
        
    #   get_vector_long()
    #-------------------------------------------------------------------
    def get_grid_long(self, var_name):

        exec("result = numpy.int32(self." + var_name + ")")
        return result
    
    #   get_grid_long()
    #-------------------------------------------------------------------
    def get_values_in_grid_long(self, var_name, IDs):

        try:
            exec("result = self." + var_name + '.flat[IDs]')
            return numpy.int32(result)
        except:
            print 'ERROR in CSDMS_base.get_values_in_grid_long().'
            print '    Returning zeros.'
            return numpy.zeros(len(IDs), dtype='int32')
        
    #   get_values_in_grid_long()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def set_scalar_long(self, var_name, scalar):

        exec("self." + var_name + " = numpy.int32(scalar)")
    
    #   set_scalar_long()
    #-------------------------------------------------------------------
    def set_vector_long(self, var_name, vector):

        #--------------------------------------------------
        # Notes: First method here should be more robust.
        #        See Notes for get_vector_double().
        #--------------------------------------------------
        exec("self." + var_name + " = numpy.array(vector, dtype='int32')")

        #-----------------------------------
        # Original method (before 2/17/10)
        #-----------------------------------
        # exec("self." + var_name + " = numpy.int32(vector)")
        
    #   set_vector_long()
    #-------------------------------------------------------------------
    def set_grid_long(self, var_name, grid):

        exec("self." + var_name + " = numpy.int32(grid)")
    
    #   set_grid_long()
    #-------------------------------------------------------------------
    def set_values_in_grid_long(self, var_name, IDs, values):

        #----------------------------------------------------------
        # Note: Type of "values" should already be long (SIDL).
        #----------------------------------------------------------
        # exec("self." + var_name + ".flat[IDs] = numpy.int32(values)")
        exec("self." + var_name + ".flat[IDs] = values")
        
    #   set_values_in_grid_long()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def get_input_items(self):

        #-------------------------------------------------------------
        # Note: There may be a way to retrieve this list auto-
        #       matically using code similar to self.has_variable().
        #-------------------------------------------------------------
        items = ['None']
        return numpy.array(items)   # (string array vs. list)
    
    #   get_input_items()
    #-------------------------------------------------------------------
    def get_output_items(self):

        items = ['None']
        return numpy.array(items)   # (string array vs. list)
    
    #   get_output_items()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def read_gui_info_file(self, gui_info_file):

        #------------------------------------------
        # Open file to read and skip header lines
        #------------------------------------------
        file_unit = open(gui_info_file, 'r')
        ## line1 = file_unit.readline()
        cfg.skip_header( file_unit, n_lines=4 )

        group_names = []
        group_sizes = []
        group_size  = 0
        #------------------
        var_names   = []
        labels      = []
        types       = []
        values      = []
        min_vals    = []
        max_vals    = []
        desc        = []

        #-----------------------------
        # Read information from file
        #-----------------------------
        while (True):
            line = file_unit.readline()
            if (line == ''):
                break
            words = line.split('|')
            char1 = words[0][0]  # (don't use .strip() due to "\n")
            if (len(words) > 6) and (char1 != '#'):
                group_size += 1  #####
                #--------------------------------------
                var_names.append( words[0].strip() )
                labels.append( words[1].strip() )
                #--------------------------------------
                vtype   = words[2].strip().lower()
                value   = words[3].strip()
                min_val = words[4].strip()
                max_val = words[5].strip()
                if (vtype != 'string'):
                    value   = eval( value )
                    min_val = eval( min_val )
                    max_val = eval( max_val )
                if (vtype in ['int', 'long']):
                    #-------------------------------
                    # Need this to avoid a warning
                    #-------------------------------
                    value   = numpy.int32( value )
                    min_val = numpy.int32( min_val )
                    max_val = numpy.int32( max_val )
                types.append( vtype )
                values.append( value )
                min_vals.append( min_val )
                max_vals.append( max_val )
                #--------------------------------------
                desc.append( words[6].strip() )
            elif (len(words) == 2) and \
                 (words[0].strip().upper() == 'PPF_GROUP_NAME'):
                group_names.append( words[1].strip() )
                if (group_size > 0):
                    group_sizes.append( group_size )
                    group_size = 0
        group_sizes.append( group_size ) # (last one)

        #--------------
        # For testing
        #--------------
##        print 'group_names =', group_names
##        print 'group_sizes =', group_sizes
##        print ' '
##        print 'var_names =', var_names
##        print 'labels    =', labels
##        print 'types     =', types
##        print 'values    =', values
##        print 'min_vals  =', min_vals
##        print 'max_vals  =', max_vals
            
        file_unit.close()
    
        return var_names, labels, types, values, min_vals, max_vals, \
               desc, group_names, group_sizes

    #   read_gui_info_file()
    #-------------------------------------------------------------------
    def get_user_input( self, services, d_services,
                        button_label="config",
                        dialog_title="Component Parameters",
                        gui_info_file=None ):

        #------------------------------------------------------
        # Note: This uses a CCA "Parameter Port" to launch a
        #       GUI dialog that prompts for input parameters
        #       to the current component instance.

        #       A call to this method can be inserted into
        #       the setServices() method of a CCA component's
        #       impl file.
        #------------------------------------------------------
        import gov.cca.ports.ParameterPortFactory
        
        #----------------------------------
        # Read lists from "gui_info_file"
        #----------------------------------
        if (gui_info_file == None):
            print 'ERROR in CSDMS_base.get_user_input():'
            print '    gui_info_file not found.'
            return
        var_names, labels, types, \
            values, min_vals, max_vals, \
            desc, group_names, group_sizes = \
            self.read_gui_info_file( gui_info_file )
        
        #----------------------------------------------
        # self.userinput = d_services.createTypeMap() 
        # Do we really need to store this in self ?
        # Seems it is only a local variable.
        # (self.userinput -> dialog)
        #----------------------------------------------
        dialog = d_services.createTypeMap()    

        try:
            port = d_services.getPort("ppf")
        except sidl.BaseException._Exception, e:
            port = None

        if (self.DEBUG):
            if (port == None):
                print 'FAILURE: Unable to get generic CCA port'
            else:
                print 'SUCCESS: Got a generic CCA port.'
        
        if not port:
            print 'getPort("ppf") returned nil'
        else: 
            ppf = gov.cca.ports.ParameterPortFactory.ParameterPortFactory( port ) 
            if not ppf:
                print 'Error casting to gov.cca.ports.ParameterPortFactory'
            else:
                if (self.DEBUG): print 'SUCCESS: Cast port to PPF port.'
            
        ppf.initParameterData(dialog, button_label)
        ppf.setBatchTitle(dialog, dialog_title)

        if (self.DEBUG):
            print 'Starting for loop over parameter list...'
            
        #----------------------------------------------
        # Add rows to prompt for requested parameters
        #--------------------------------------------------------
        # Note that var_names must be unique or they won't show
        #--------------------------------------------------------
        group_num = 0
        group_sum = 0
        for k in xrange(len(var_names)):
##            if (self.DEBUG):
##                print 'k =', k   #############
##                print 'types[k].lower() =', types[k].lower()
            
            #----------------------------------------------
            # Start a new "panel" with new "group name" ?
            #----------------------------------------------
            if (k == group_sum):
                ppf.setGroupName( dialog, group_names[ group_num ] )
                group_sum += group_sizes[ group_num ]
                group_num += 1
##                print 'group_num =', group_num
##                print 'group_sum =', group_sum
                
            #------------------------------------------
            # Prompt for a variable of requested type
            #------------------------------------------
            var_name = var_names[k]
            var_type = types[k].lower()
            if   (var_type in ['double', 'float']):
                ppf.addRequestDouble(dialog, var_names[k],
                                     desc[k], labels[k], values[k],
                                     min_vals[k], max_vals[k])
            elif (var_type in ['int', 'long']):
                ppf.addRequestInt(dialog, var_names[k],
                                  desc[k], labels[k], values[k],
                                  min_vals[k], max_vals[k])
            elif (var_type == 'string'):
                ppf.addRequestString(dialog, var_names[k],
                                     desc[k], labels[k], values[k] )
                                     ### min_vals[k], max_vals[k])
            else:
                print 'ERROR in CSDMS_base.get_user_input().'
                print '   Unsupported type name: ' + '__' + types[k] + '__'
                
        if (self.DEBUG):
            print 'Finished with for loop.'
            print 'Calling ppf.addParameterPort()...'
        
        ppf.addParameterPort(dialog, services)
        d_services.releasePort("ppf")  
        if (self.DEBUG):
            #-------------------------------------------------
            # Note: This message is printed even if the user
            #       has not yet opened the PPF dialog.
            #-------------------------------------------------
            print 'Finished with get_user_input().'
        
        #-----------------------------------------------
        # Does call to "releasePort()" imply a closing
        # of the dialog ???
        #-----------------------------------------------
        # Maybe need to sleep for a bit here.
        # time.sleep(0.2)
        
        #---------------------------------------------------------
        # Try to store the user-entered values directly into
        # the current component's state.  This probably requires
        # that user has dismissed dialog with "OK" button.
        #---------------------------------------------------------
        # If we can't do this here, then wrap this block into a
        # small function that can be called by "go()" method.
        #---------------------------------------------------------
        self.dialog = dialog
        ## self.load_user_input(var_names, types, dialog_str='dialog')

        #--------------------------------------------------------
        # Store these for use by "CSDMS_base.load_user_input()"
        #--------------------------------------------------------
        # Note: Last line doesn't work quite as intended because
        #       this function is called (and DEBUG message is
        #       printed, see above) even if user has not yet
        #       opened the PPF dialog. (But need it for now.)
        #----------------------------------------------------------  
        self.PPF_var_names = var_names
        self.PPF_types     = types
        self.USER_SET_VALUES = True
        
    #   get_user_input()
    #-------------------------------------------------------------------
    def load_user_input(self, var_names=None, types=None,
                        dialog_str='self.dialog'):

        if (self.DEBUG):
            print 'Saving user input into component state.'

        if (var_names == None):
            var_names = self.PPF_var_names
        if (types == None):
            types = self.PPF_types

        #-----------------------------------------
        # Save user input into component's state
        #-----------------------------------------
        for k in xrange(len(var_names)):
            #-----------------------------------------------
            # Get the variable's data type.  Allowed types
            # are:  int, long, float, double, string, file.
            # "file" is special and should not be used in
            # the GUI_info_file (see below).
            #-----------------------------------------------
            var_name  = var_names[k]
            var_type  = types[k].lower()
            
            if (self.has_variable( var_name + "_type" )):  
                #-----------------------------------------
                # This requires that the variable's type
                # was requested in a preceding row and
                # has already been saved into self.
                #-----------------------------------------
                # TF types are Scalar, Time_Series, Grid
                # and Grid_Sequence.
                #-----------------------------------------
                exec("TF_type = self." + var_name + "_type")
                if (TF_type.lower() == 'scalar'):
                    var_type = 'double'
                else:
                    var_type = 'file'
                
            #------------------------------------
            # Read variable with that data type
            #------------------------------------            
            if (var_type in ['double', 'float']):            
                exec( "self." + var_name + " = " +
                      dialog_str + ".getDouble( var_name, 0.0 )" )
            elif (var_type in ['long', 'int']):                
                exec( "self." + var_name + " = " +
                      dialog_str + ".getInt( var_name, 0 )" )
            elif (var_type == 'string'):
                exec( "self." + var_name + " = " +
                      dialog_str + ".getString( var_name, '0.0' )" )
            elif (var_type == 'file'):
                #-------------------------------------------
                # Append "_file" to var_name, then save it
                #------------------------------------------- 
                exec( "self." + var_name + "_file = " +
                      dialog_str + ".getString( var_name, '0.0' )" ) 
            else:
                print 'ERROR in CSDMS_base.load_user_parameters().'
                print '   Unsupported data type = ' + types[k] + '.'
                print ' '
    
    #   load_user_input()
    #-------------------------------------------------------------------
    def go(self):

        self.run_model()

    #   go()
    #-------------------------------------------------------------------
    def run_model(self, directory=None, data_prefix=None,
                  case_prefix=None, mode="main",
                  n_steps=5):

        #------------------------------------------------------
        # If directory, etc. were set with the GUI, then use
        # those settings.  Otherwise, use functions shown.
        #------------------------------------------------------
        if (directory == None):
            if (self.directory != None):
                directory = self.directory
            else:
                directory   = tf_utils.TF_Test_Directory()
                # directory   = tf_utils.Current_Directory()
        #------------------------------------------------------
        if (data_prefix == None):
            if (self.data_prefix != None):
                data_prefix = self.data_prefix
            else:
                data_prefix = tf_utils.TF_Test_Data_Prefix()
        #------------------------------------------------------
        if (case_prefix == None):
            if (self.case_prefix != None):
                case_prefix = self.case_prefix
            else:
                case_prefix = tf_utils.TF_Test_Case_Prefix()
                # case_prefix = tf_tuils.Get_Case_Prefix()

        #---------------------------
        # Initialize the model run
        #---------------------------
        self.initialize(directory=directory,
                        case_prefix=case_prefix,
                        data_prefix=data_prefix, mode=mode)###

        #---------------------------------------------------------- 
        # Note:  If a component calls "initialize_time_vars()" in
        #        its "initialize()" method, then "self.DONE" will
        #        be set to False.
        #----------------------------------------------------------
        #        Should it be set to False in "__init__()" just
        #        in case "initialize_time_vars()" is not called
        #        or is over-ridden??
        #----------------------------------------------------------
        if not(hasattr(self, 'DONE')):
            self.DONE = False
            
        #----------------------------------------------------------- 
        # Note:  If the number of timesteps is specified in a
        #        component's CFG file and is then saved by
        #        "read_config_file()" as "n_steps" then we should
        #        honor that setting.  Otherwise we use the n_steps
        #        argument.
        #-----------------------------------------------------------    
        if (hasattr(self, 'n_steps')):
            ## print 'NUMBER OF STEPS =', self.n_steps  ####
            n_steps = self.n_steps
                
        while not(self.DONE):
            if (self.DEBUG):
                #-------------------------------------------
                # Exceptions will not be caught and
                # Python error messages will be displayed.
                #-------------------------------------------
                print 'time_index =', self.time_index
                self.update()
            else:   
                try:
                    self.update()
                except:
                    print 'ERROR in run_model() method at:'
                    print '   time_index =', self.time_index
                    self.status = 'failed'
                    self.DONE = True

            #------------------------------------------------
            # If the model has set self.DONE = True, then
            # stop, even if we haven't reached n_steps yet.
            #------------------------------------------------
            TIMES_UP  = (self.time_index >= n_steps)
            self.DONE = (self.DONE or TIMES_UP)

        #-------------------------
        # Finalize the model run
        #-------------------------
        self.finalize()

    #   run_model()
    #-------------------------------------------------------------------
##    def read_config_file(self):
##
##        #------------------------------------------
##        # Read parameters from a CFG file that is
##        # in the current working directory.
##        #------------------------------------------
##        print 'Process component: Reading config file...'
##        self.config_file = (self.directory +
##                            self.case_prefix + '_PROC.cfg')
##        file_unit = open(self.config_file, 'r')
##
##        #-----------------------------
##        # Skip over the header lines
##        #-----------------------------
##        for k in xrange(1, 5):
##            line = file_unit.readline()
##            
##        #------------------------
##        # Read the process vars
##        #------------------------
####        method       = Read_Vars(file_unit, data_type='BYTE')
####        method_name  = Read_Vars(file_unit, data_type='STRING')
####        n_layers     = Read_Vars(file_unit, data_type='INTEGER')
####        dt_type, dt  = Read_Vars(file_unit, True)
####        #-----------------------------------------
####        self.method      = method
####        self.method_name = method_name
####        self.n_layers    = n_layers
####        self.dt          = float64(dt)   # [seconds]
## 
##        #-----------------------
##        # Close the config file
##        #-----------------------
##        file_unit.close()
##
##        #---------------------------------------------------------
##        # Make sure that all "save_dts" are larger or equal to
##        # the specified process dt.  There is no point in saving
##        # results more often than they change.
##        # Issue a message to this effect if any are smaller ??
##        #---------------------------------------------------------
####        self.save_grid_dt    = maximum(self.save_grid_dt,    self.dt)
####        self.save_pixels_dt  = maximum(self.save_pixels_dt,  self.dt)
##        
##    #   read_config_file()
    #-------------------------------------------------------------
    def initialize(self, directory=None,
                   data_prefix=None,
                   case_prefix=None, mode="module"):

        self.status = 'initializing'
        self.set_directory(directory, data_prefix, case_prefix)
        self.read_grid_info()
       
        #---------------------------------------------------
        # Load component parameters from a config file ?
        #---------------------------------------------------      
##        self.read_config_file()
##        self.get_ports()
        self.initialize_time_vars()

        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
##        self.open_input_files()
##        self.read_input_files()

        #-------------------------------
        # Save these for writing output
        #-------------------------------
##        self.outlet_ID  = self.bp.get_scalar_long('outlet_ID')
##        n_outlets       = self.bp.get_size('outlet_IDs')
##        self.outlet_IDs = self.bp.get_ID_list('outlet_IDs', n_outlets)
##        self.open_output_files()

        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
##    def update(self, time_seconds=None):
##
##        self.status = 'updating'  # (OpenMI 2.0 convention)
##
##        #-----------------------------------
##        # Update computed values here with
##        # a series of "update_var()" calls
##        #----------------------------------
##        # self.update_var1()
##        # self.update_var2()
##
##        #------------------------
##        # Update internal clock
##        #------------------------
##        self.update_time()
##
##        #-------------------------------
##        # Check for NaNs, etc. in var1
##        #-------------------------------    
##        # self.check_var1()
##
##        #------------------------------------------
##        # Read next infil vars from input files ?
##        #------------------------------------------
##        self.read_input_files()
##
##        #----------------------------------------------
##        # Write user-specified data to output files ?
##        #----------------------------------------------
##        self.write_output_files(time_seconds)
##        self.status = 'updated'  # (OpenMI 2.0 convention)
##        
##    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI 2.0 convention)
        self.close_input_files()    ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI 2.0 convention)
        
    #   finalize()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def initialize_required_components(self, mode="module"):

##        print "MODE     =", mode
##        print "self.CCA =", self.CCA
##        print ' '
        
        if (mode == "main"):
            if not(self.CCA):
                #-----------------------------------------
                # Create embedded "process object ports"
                #-----------------------------------------
                self.embed_child_components()
                self.add_child_ports()
            if not(self.SILENT):
                print 'INITIALIZING CCA PORTS....'
                print ' '
            self.initialize_ports()
        elif (mode == "module"):
            #-----------------------------------------------
            # In "initialize()" method of component's Impl
            # file it must call "get_cca_ports() before
            # calling "initialize()".
            #-----------------------------------------------
            pass
        
    #   initialize_required_components()
    #-------------------------------------------------------------------
    def get_cca_port(self, port_name, port_type, project_name,
                     d_services):

        exec("import " + project_name)
        
        #--------------------------------------------
        # Note:  Use "d_services" added by Bocca to
        #        an Impl file to get a CCA port
        #--------------------------------------------
        SUCCESS = True
        str2 = project_name + "." + port_type + "." + port_type
        
        #-----------------------------------
        # Try to get a "generic" CCA port.
        #-----------------------------------
        try:
            port = d_services.getPort( port_name )
            print "SUCCESS: Got CCA port: " + port_name
            ## print "*** type(port) =", type(port)
        except:
            port    = None
            message = "FAILURE:  Unable to get CCA port: " + port_name
            print message
            print >> sys.stderr, message
            SUCCESS = False
            
        #------------------------------------------
        # Try to typecast the port to "port_type"
        #------------------------------------------
        exec( "new_port = " + str2 + "( port )") 

        if (new_port == None):
            print 'FAILURE: Unable to cast CCA port: ' + port_name
            d_services.releasePort( port_name )
            SUCCESS = False

        return new_port

    #   get_cca_port()  
    #-------------------------------------------------------------------
    def get_cca_ports(self, port_names, short_names,
                      port_type, project_name, d_services):

        exec("import " + project_name)
        
        #--------------------------------------------
        # Note:  Use "d_services" added by Bocca to
        #        an Impl file to get a CCA port
        #--------------------------------------------
        SUCCESS = True
        str2 = project_name + "." + port_type + "." + port_type
        str3 = "( port )"
        
        for k in xrange(len(port_names)):
            #-----------------------------------
            # Try to get a "generic" CCA port.
            #-----------------------------------
            name = port_names[k]
            try:
                port = d_services.getPort( name )
                print "SUCCESS: Got CCA port: " + name
                ## print "*** type(port) =", type(port)
            except:
                port    = None
                message = "FAILURE:  Unable to get CCA port: " + port_name
                print message
                print >> sys.stderr, message
                SUCCESS = False
                
            #------------------------------------------
            # Try to typecast the port to "port_type"
            # and then store it within "self"
            #------------------------------------------
            str1 = "self." + short_names[k]  
            exec( str1 + " = " + str2 + str3 )
 
            exec( "UNABLE = (" + str1 + " == None)" )
            if (UNABLE):
                print 'FAILURE: Unable to cast CCA port: ' + name
                exec( "d_services.releasePort( "  + name + " )")
                SUCCESS = False

        return SUCCESS
    
        #-------------------------------------------------------------
        # An example where:
        #    port_name    = "channels"
        #    short_name   = "cp"
        #    port_type    = "IRFPort"
        #    project_name = "topoflow"

        #    In CCA Impl file before call to initialize():
        
        #    OK = self.tf.add_cca_ports( ["channels"], ["cp"],
        #                                "IRFPort", "topoflow",
        #                                 self.d_services )
        #    if not(OK): return -1
        #-------------------------------------------------------------
        # port = self.d_services.getPort( "channels" )
        # self.cp = topoflow.IRFPort.IRFPort( port )        
        # if (self.cp == None):
        #     print "Surface flow component is not connected."
        #     self.d_services.releasePort( "channels" )
        #     return -1
        #-------------------------------------------------------------
        # Example call in initialize method of "TopoFlow_Impl.py":
        #
        # port_names = ['channels', 'precip', 'snow', 'evap',
        #               'infil', 'sat_zone', 'diversions',
        #               'meteorology', 'basins']
        # short_names = ['cp', 'pp', 'sp', 'ep', 'ip', 'gp',
        #                'dp', 'mp', 'bp']
        # OK = self.tf.get_cca_ports( port_names, short_names,
        #                             "IRFPort", "topoflow",
        #                             self.d_services )
        # if not(OK): return -1
        #-------------------------------------------------------------
        
    #   get_cca_ports()
    #-------------------------------------------------------------------
    def release_cca_ports(self, port_names, d_services):
   
        #----------------------
        #  Release all ports
        #----------------------
        for name in port_names:
            d_services.releasePort( name )
                 
    #   release_cca_ports()    
    #-------------------------------------------------------------------
    def add_child_port(self, port1_str, port2_str, SELF=False):

        #-----------------------------------
        # Example:  port1_str = "cp"
        #           port2_str = "pp"
        #     =>    self.cp.pp = self.pp
        #--------------------------------------
        # Note: If (self == channels (cp)),
        #       then want "self.ep.cp = self"
        #--------------------------------------
        str1 = "self." + port1_str + "."
        str2 = port2_str + " = "
        if not(SELF):
            str3 = "self." + port2_str
        else:
            str3 = "self"
        exec( str1 + str2 + str3)
            
    #   add_child_port()
    #-------------------------------------------------------------------
    def get_port_data(self, var_name, port=None,
                      port_name='CCA', VERBOSE=False):
        
        #---------------------------------------------------------
        # Note: This method is "private", or not part of the
        #       exposed component interface (or port), so its
        #       return type (in Python) can be dynamic.
        #       However, the functions called by this one below
        #       are "port functions" and therefore have static
        #       return types.
        #---------------------------------------------------------
        # Note: To maximize performance, we may want to just get
        #       the type for each var_name once during the
        #       initialize() call and save it internally.
        #---------------------------------------------------------
        if (port is None):
            print "ERROR: get_port_data() needs CCA port argument."
            return numpy.float64(0)
        
        #---------------------------------------------------------
        try:
            if (port.is_scalar(var_name)):
                data = port.get_scalar_double( var_name )
                str2 = 'SCALAR'
            elif (port.is_grid(var_name)):
                data = port.get_grid_double( var_name )
                str2 = 'GRID'
            elif (port.is_vector(var_name)):
                data = port.get_vector_double( var_name )
                str2 = 'VECTOR'
            else:
                print 'ERROR in CSDMS_base.get_port_data().'
                print '    ' + var_name + ' is not SCALAR, VECTOR or GRID.'
                print '    Returning scalar value of 0.'
                print ' '
                data = numpy.float64(0)
                str2 = 'UNKNOWN_TYPE'
        except:
            data = numpy.float64(0)
            str2 = 'SCALAR'
            print 'ERROR in CSDMS_base.get_port_data().'
            pn_str = ' from ' + port_name + ' port.'
            print '   Could not get ' + var_name + pn_str
            print '   Returning scalar value of 0.'
            print ' '
             
        #--------------
        # For testing
        #--------------
        if (VERBOSE):
            str1 = '>>> ' + var_name + ' is a '
            str3 = '. type(' + var_name + ') =', type(data)
            print (str1 + str2 + str3)
             
        #----------------------------------
        # Make sure return type is double
        #----------------------------------
        return numpy.float64(data)
        
    #   get_port_data()
    #-------------------------------------------------------------------
    def get_port_data_OLD(self, var_name, port=None):
        
        #---------------------------------------------------------
        # Note: This method is "private", or not part of the
        #       exposed component interface (or port), so its
        #       return type (in Python) can be dynamic.
        #       However, the functions called by this one below
        #       are "port functions" and therefore have static
        #       return types.
        #---------------------------------------------------------
        # Note: To maximize performance, we may want to just get
        #       the type for each var_name once during the
        #       initialize() call and save it internally.
        #---------------------------------------------------------
        if (port is None):
            print "ERROR: get_port_data() needs CCA port argument."
            return -1
        #---------------------------------------------------------        
##        if (port is None):
##            exec("return self." + var_name)
##        else:
        #---------------------------------------------------------
        if (port.is_scalar(var_name)):
            # print '>>> ' + var_name + ' is a SCALAR.'
            ## return port.get_scalar_double( var_name )
            scalar = port.get_scalar_double( var_name )
            # print '>>> type(' + var_name + ') =', type(scalar)
            return numpy.float64( scalar )
        elif (port.is_grid(var_name)):
            # print '>>> ' + var_name + ' is a GRID.'
            ## return port.get_grid_double( var_name )
            grid = port.get_grid_double( var_name )
            # print '>>> type(' + var_name + ') =', type(grid)
            return numpy.float64( grid )
        else:
##            nv = port.get_size( var_name )
##            if (nv == 1):
##                return numpy.float64( 
            print 'ERROR in CSDMS_base.get_port_data().'
            print '    Variable is not SCALAR or GRID.'
            print '    rank(var) =', port.get_rank( var_name )
            print '    Returning 0.'
            print ' '
            return numpy.float64(0)
            
    #   get_port_data_OLD()
    #-------------------------------------------------------------------
    def set_port_data(self, var_name, var, port=None):
        
        #---------------------------------------------------------
        # Note: This method is "private", or not part of the
        #       exposed component interface (or port), so its
        #       return type (in Python) can be dynamic.
        #       However, the functions called by this one below
        #       are "port functions" and therefore have static
        #       return types.
        #---------------------------------------------------------
        # Note: To maximize performance, we may want to just get
        #       the type for each var_name once during the
        #       initialize() call and save it internally.
        #---------------------------------------------------------
        # Note: (2/8/10) Bug fix.  A variable such as "h_swe"
        #       may be initialized to a scalar (in snow comp.)
        #       so "port.is_scalar()" will return True.  But if
        #       "var" is a grid, then we won't be allowed to
        #       "set it" as one.  Solution is to replace:
        #       if (port.is_scalar(var_name)) to:
        #       if (size(var) == 1):
        #---------------------------------------------------------            
##        if (port is None):
##            exec("return self." + var_name)
##        else:
        #----------------------------------------------
        if (size(var) == 1):
            port.set_scalar_double(var_name, var)
        else:
            port.set_grid_double(var_name, var)
        #----------------------------------------------
##        if port.is_scalar(var_name):
##            port.set_scalar_double(var_name, var)
##        else:
##            port.set_grid_double(var_name, var)

    #   set_port_data()
    #-------------------------------------------------------------------
    def get_rank(self, var_name):

        exec("rank = numpy.rank(self." + var_name + ")")
        return rank
    
    #   get_rank()
    #-------------------------------------------------------------------
    def get_size(self, var_name):

        #-------------------------------------------------------
        # Notes: This is used by a caller to determine the
        #        number of elements a given variable has.
        #        This information can then be used to size
        #        an array, for example.  See "get_rank()".

        #        In a dynamically-typed language like Python,
        #        the dynamic typing can be used with NumPy to
        #        allow very flexible input types.

        # NB!    Right now, case in var_name must be an exact
        #        match.
        #-------------------------------------------------------
        exec("n = numpy.size(self." + var_name + ")")
        return n
    
    #   get_size()
    #-------------------------------------------------------------------
    def set_directory(self, directory, data_prefix, case_prefix):

        #-------------------------------------
        # Make sure these are always defined
        #-------------------------------------
        if (directory is None):
            directory = tf_utils.Current_Directory()
        if (data_prefix is None):
            data_prefix = 'Null'
            ## Get_Data_Prefix has filepath arg and
            ## a different purpose.
            ## data_prefix = tf_utils.Get_Data_Prefix()
        if (case_prefix is None):
            case_prefix = tf_utils.Get_Case_Prefix()

        #--------------------------------------------------
        # Add trailing separator to directory, if missing
        #--------------------------------------------------
        if (directory[-1] != os.sep):
            directory += os.sep
        
        self.directory   = directory
        self.data_prefix = data_prefix
        self.case_prefix = case_prefix

        #----------------------------
        # CD to working directory ?
        # May contain blank spaces
        #----------------------------
        os.chdir( self.directory )
        
    #   set_directory()
    #-------------------------------------------------------------------
    def read_grid_info(self):

        #------------------------------------------
        # Read grid info from an RTI file that is
        # in the current working directory.
        #------------------------------------------
        if not(self.SILENT):
            print 'Process component: Reading grid info...'
        self.grid_info_file = (self.directory +
                               self.data_prefix + '.rti')


        info = rti_files.read_info( self.grid_info_file )
        
        #----------------------
        # Convenient synonyms
        #-----------------------------------------------------
        # Note that "info" has additional derived attributes
        # such as: n_pixels, bpe, grid_size and SWAP_ENDIAN.
        #-----------------------------------------------------
        self.rti = info
        self.nx  = info.ncols
        self.ny  = info.nrows
        
        #------------------------------------------------
        # Get grid cell areas, "da", which is either a
        # scalar (if same for all grid cells) or a grid
        # with default units of "m^2".
        #------------------------------------------------
        self.da = pixels.get_da( info )
        
##        rti = rti_file.rti_file(self.grid_info_file)
##        OK = rti.read_file()
##        if not(OK):
##            print 'ERROR: Could not read grid info from file:'
##            print self.grid_info_file
##            print ' '
##            return
##        self.rti      = rti
##        self.nx       = rti.ncols
##        self.ny       = rti.nrows
##        self.n_pixels = rti.n_pixels
##        self.da = pixels.get_da(rti)
 
    #   read_grid_info()
    #-------------------------------------------------------------------
    def store_outlet_IDs(self):
        
        outlet_IDs = self.bp.get_vector_long('outlet_IDs')
        outlet_ID  = outlet_IDs[0]
##        self.outlet_IDs = outlet_IDs   # (long-int calendar indices)
##        self.outlet_ID  = outlet_ID
        self.outlet_IDs = (outlet_IDs / self.nx, outlet_IDs % self.nx)
        self.outlet_ID  = (outlet_ID  / self.nx, outlet_ID  % self.nx)

    #   store_outlet_IDs()   
    #-------------------------------------------------------------------
    def initialize_time_vars(self, units='seconds'):

        #------------------
        # Start the clock
        #------------------
        self.start_time = time.time()
        
        #--------------------------------
        # Initialize the time variables
        #--------------------------------
        self.time_units = units.lower()
        self.time_index = int32(0)
        self.time       = float64(0)
        self.DONE       = False
        
        #--------------------------
        # Time conversion factors
        #--------------------------
        self.sec_per_year = float64(365) * 24 * 3600
        self.min_per_year = float64(365) * 24 * 60
        
        #-------------------------------------------
        # For backward compatibility with TopoFlow
        #-------------------------------------------
        self.time_sec = float64(0)
        self.time_min = float64(0)
            
        #--------------------------------------------
        # For print_time_and_value() function below
        #--------------------------------------------
        self.last_print_time = time.time()
        
##        self.last_check_time  = time.time()  # (for user interrupt)
##        self.last_plot_time   = float32(0.0)   ### CHECK ###
        
    #   initialize_time_vars()
    #-------------------------------------------------------------------
    def update_time(self):

        #---------------------
        # Increment the time
        #---------------------
        self.time_index += 1
        self.time       += self.dt  # (use same units as dt)
        
        if (self.time_units == 'seconds'):
            self.time_sec = self.time                    # [seconds]
            self.time_min = self.time_sec / float64(60)  # [minutes]
        elif (self.time_units == 'years'):
            #-----------------------------------
            # Used by GC2D and Erode (12/4/09)
            #-----------------------------------
            self.time_sec = self.time * self.sec_per_year  ####
            self.time_min = self.time_sec / float64(60)  # [minutes]
            
    #   update_time()
    #-------------------------------------------------------------------
    def print_time_and_value(self, var, var_name='Q_out',
                             units_name='[m^3/s]',
                             interval=2.0,
                             PRINT_INDEX=False):

        #---------------------------------------------------
        # Note: Print the model time, in minutes, and the
        #       current value of "var", at the specified
        #       real-time "interval" (in seconds).
        #---------------------------------------------------
        # Note: Plotting hydrograph at same interval is
        #       generally too infrequent.
        #---------------------------------------------------
        #  self.Tstr is set in TF by initialize_stop_vars()
        #---------------------------------------------------
        elapsed_time = (time.time() - self.last_print_time)
        if (elapsed_time > interval):
            if (self.time_units == 'seconds'):
                cur_time = self.time_min
                time_units_str = ' [min]'
            else:
                cur_time = self.time
                time_units_str = ' [' + self.time_units + ']' 
            time_str = 'Time = ' + ("%10.2f" % cur_time)
            time_str = time_str + time_units_str
            #-------------------------------------------------
            var_str  = var_name + ' = ' + ("%10.5f" % var)
            var_str  = var_str  + ' ' + units_name          
            #-------------------------------------------------      
            print (time_str + ',  ' + var_str)
            #-----------------------------------------------------
            if (PRINT_INDEX):
                print 'n =', self.time_index, 'of', self.n_steps
            #-----------------------------------------------------                
            self.last_print_time = time.time()

    #   print_time_and_value()    
    #-------------------------------------------------------------------
    def print_run_time(self, proc_name='component',
                       sec_digits=4, seconds=None,
                       SILENT=None):
                       ### SUB_PROCESS=False) 

        #------------------------------------------------------
        # If "seconds" argument is only provided for testing.
        # You can provide this value to make sure that the
        # minuts, hours, days, etc. are computed correctly.
        #------------------------------------------------------
        if (seconds == None):    
            finish  = time.time()
            seconds = (finish - self.start_time)

        #----------------------------------
        # Compute minutes, hours and days
        #----------------------------------
        dec_part  = (seconds % float32(1.0))     #(Save decimal part)
        days      = int32(seconds) / int32(86400)
        secs_left = int32(seconds) % int32(86400)
        hours     = (secs_left / int32(3600))
        secs_left = (secs_left % int32(3600))
        minutes   = (secs_left / int32(60))
        seconds   = (secs_left % int32(60))
        #-----------------------------------------
        #hours     = long(seconds)  /  3600L
        #secs_left = long(seconds) mod 3600L
        #minutes   = (secs_left  /  60L)
        #seconds   = (secs_left mod 60L)
        
        #----------------------------
        # Construct the time string
        #----------------------------
        time_string = ''
        #--------------------------------------------------------
        if (days > 0):    
            if (days > 1):    
                e0 = ' days, '
            else:    
                e0 = ' day, '
            time_string += str(days) + e0
        #--------------------------------------------------------
        if (hours > 0):    
            if (hours > 1):    
                e1 = ' hours, '
            else:    
                e1 = ' hour, '
            time_string += str(hours) + e1
        #--------------------------------------------------------
        if (minutes > 0):    
            if (minutes > 1):    
                e2 = ' minutes, '
            else:    
                e2 = ' minute, '
            time_string += str(minutes) + e2
        
        #-----------------------------------------
        # Default is 4 digits after the decimal.
        #-----------------------------------------
        dec_pastr = ('.' + str(dec_part)[2:2+sec_digits])
        time_string += str(seconds) + dec_pastr + ' seconds.'
        
        if not(SILENT):
            print ('Run time for ' + proc_name + ' = ')
            print  time_string
            print ' '
            
##            if (SUB_PROCESS):    
##                PART1 = '>> '
##            else:    
##                PART1 = ''
##            print (PART1 + 'Run time for ' + procname + ' = ')
##            print (PART1 + time_string)
##            print ' '
     
    #   print_run_time()
    #-------------------------------------------------------------------    
