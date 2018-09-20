'''This file defines inputs by pulling nk from files'''



from numpy import  inf, loadtxt, pi

import re






class dyn_basic_setup():
    def __init__(self,thickness=10,R_line='0 50 60',R_dir='Reflection_10nm_CuAg_on_silica_substrate.txt',
                    T_dir='Transmission_10nm_CuAg_on_silica_substrate.txt',T_line=None,):
        
    #for identity
    #thickness=self.le_t.text()
    
      
        self.layer_index_of_fit=1
    
        self.lamda_min = 300
        self.lamda_max = 1000
    
        kind = 'cubic'
    
        #### air layer
      
        nk_f_list = [self.nk_f_air] #top working down
       
        self.thickness_list = [inf]
        coherency_list = ['i']
        
        ##### film layers  ######
        film_thickness = thickness # nm
        #nk_f_list.append(functionize_nk_file('Au-glass_10nm_30p_effective_nk.txt',skiprows=1, kind = kind))
        nk_f_list.append(self.nk_f_air)
        self.thickness_list.append(film_thickness)
        coherency_list.append('c')
        
        ###### substrate layer
        substrate_thickness = 0.5e-3 *(1e9) # 2 mm in nm
        nk_f_list.append(self.nk_f_silica)
        self.thickness_list.append(substrate_thickness)
        coherency_list.append('i')
        
        ##### back air layer ####
        nk_f_list.append(self.nk_f_air)
        self.thickness_list.append(inf)
        coherency_list.append('i')
        
        ###########
        fit_nk_f = nk_f_list[self.layer_index_of_fit]
        
        
        #################  experimental data
        #fetch input from LineEdit
        R=re.split(r"[\s,]+",R_line)
    
        
      
        self.spectrum_function_list = []
        spectrum_name_list = [] # this is for me to label things later
        
        with open(R_dir,'r') as f:
            R_index=f.readline()
        R_index=re.split(r" *\t+ *| {2,}",R_index)
       
        R_data = loadtxt(R_dir, skiprows = 2).T
        lamda = R_data[0]
        self.lamda_min=lamda[0]
        self.lamda_max=lamda[len(lamda)-1]
        index=3   #in case not contain the input degree
        for R_value in R:
            if R_value=='0':
                self.spectrum_function_list.append( self.extrap(lamda, R_data[1]/100.0, kind = 'linear' ) )
                spectrum_name_list.append('0 deg Reflection')
            else:
                pattern=re.compile(r".*"+R_value+".*")       #use re.compile to add variable in a RegE
                for i in range(len(R_index)):
                    if re.match(pattern, R_index[i])!=None:
                        index=i
                        break
       
                self.spectrum_function_list.append( self.extrap(lamda, R_data[(index-2)*2+2]/100.0, kind = 'linear' ) )
                spectrum_name_list.append('%s deg S-polarization Reflection'%R_value)
                self.spectrum_function_list.append( self.extrap(lamda, R_data[(index-2)*2+3]/100.0, kind = 'linear' ) )
                spectrum_name_list.append('%s deg P-polarization Reflection'%R_value)
            
         # lnear interpoation prevents interpoation of TRA values outside 0-100%
        
        T_data = loadtxt(T_dir, skiprows = 1).T
        self.spectrum_function_list.append(  self.extrap(T_data[0], T_data[1]/100.0, kind = 'linear' ) )
        spectrum_name_list.append('0 deg Transmission')
       
       
       
       
        #define the parameter list
      
        self.thickness_list[self.layer_index_of_fit] = thickness
        lamda=300
        # order must match the spectrum_list_generator
       
        self.param_list = []
        for R_value in R:
            self.param_list.append( {
                'lamda' : lamda,
                'snell_angle_front' : float(R_value) * pi/180,
                'layer_index_of_fit' : self.layer_index_of_fit,
                'nk_f_list' : nk_f_list,
                'thickness_list' : self.thickness_list,
                'coherency_list' : coherency_list,
                'tm_polarization_fraction' : 0.0,
                'spectrum' : 'R'} )
            if R_value!='0':
                self.param_list.append( {
                'lamda' : lamda,
                'snell_angle_front' : float(R_value) * pi/180,
                'layer_index_of_fit' : self.layer_index_of_fit,
                'nk_f_list' : nk_f_list,
                'thickness_list' : self.thickness_list,
                'coherency_list' : coherency_list,
                'tm_polarization_fraction' : 1.0,
                'spectrum' : 'R'} )
    
    
        self.param_list.append( {
                'lamda' : lamda,
                'snell_angle_front' : 0 * pi/180,
                'layer_index_of_fit' : self.layer_index_of_fit,
                'nk_f_list' : nk_f_list,
                'thickness_list' : self.thickness_list,
                'coherency_list' : coherency_list,
                'tm_polarization_fraction' : 0.0,
                'spectrum' : 'T'} )

    
    def nk_f_air(self,lamda):
        return 1.0+0.0j*lamda
            
    def nk_f_silica(self,lamda):
        return 1.5+0.0j*lamda
  
    ### this function is what we are creating as
    def spectrum_list_generator(self,lamda): # for measured or simulated spectra
    
        spectrum_list = [spectrum_function(lamda) for spectrum_function in self.spectrum_function_list] # why not just directly use the spectrum_function_list?
        #because this allows us in the future to use spectra with different Wavelength ranges!
        #return spectrum_list #must be a list of spectra  matching the params
    
        return spectrum_list
        
        
        #################################### illumination conditions
    def parameter_list_generator(self,lamda,thickness=None): # layer geometry and neighbor properties
    
        if thickness!=None:
            self.thickness_list[self.layer_index_of_fit]=thickness
            for param in self.param_list:
                param['thickness_list']=self.thickness_list
        for param in self.param_list:
                param['lamda']=lamda
        
        return self.param_list
    
    def get_lamda_range(self):
        return self.lamda_min,self.lamda_max
    def extrap(self,lamda, n, kind = 'linear'):
        '''Requires that lamda be in increasing order'''
        upper_value = n[-1]
        lower_value = n[0]
        from scipy.interpolate import interp1d, BSpline
        def is_in_bounds(self,lamda):
            if (self.lower_bound <= lamda) and (lamda <= self.upper_bound):
                return True
            else:
                return False
    
        # now we instantiate
        if kind != 'cubic_bspline':
            interp1d.upper_bound = 0
            interp1d.lower_bound = 0
            interp1d.is_in_bounds = is_in_bounds
            func = interp1d(lamda, n, kind=kind, bounds_error = False, fill_value = (lower_value, upper_value))
        else:
            BSpline.upper_bound = 0
            BSpline.lower_bound = 0
            BSpline.is_in_bounds = is_in_bounds
            func = BSpline(lamda, n, k = 3, extrapolate = True)
    
        func.upper_bound = max(lamda)
        func.lower_bound = min(lamda)
        return func  
