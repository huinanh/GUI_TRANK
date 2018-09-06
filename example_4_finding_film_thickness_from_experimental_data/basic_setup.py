
'''This file defines inputs by pulling nk from files'''


from TRANK import functionize_nk_file, extrap
from numpy import  inf, arange, loadtxt, pi



layer_index_of_fit = 1

lamda_min = 300
lamda_max = 1000

kind = 'cubic'

#### air layer
def nk_f_air(lamda):
	return 1.0+0.0j*lamda

def nk_f_silica(lamda):
	return 1.5+0.0j*lamda



nk_f_list = [nk_f_air] #top working down
thickness_list = [inf]
coherency_list = ['i']

##### film layers  ######
film_thickness = 10 # nm
#nk_f_list.append(functionize_nk_file('Au-glass_10nm_30p_effective_nk.txt',skiprows=1, kind = kind))
nk_f_list.append(nk_f_air)
thickness_list.append(film_thickness)
coherency_list.append('c')

###### substrate layer
substrate_thickness = 0.5e-3 *(1e9) # 2 mm in nm
nk_f_list.append(nk_f_silica)
thickness_list.append(substrate_thickness)
coherency_list.append('i')

##### back air layer ####
nk_f_list.append(nk_f_air)
thickness_list.append(inf)
coherency_list.append('i')

###########
fit_nk_f = nk_f_list[layer_index_of_fit]


#################  experimental data

spectrum_function_list = []
spectrum_name_list = [] # this is for me to label things later

R_data = loadtxt('Reflection_10nm_CuAg_on_silica_substrate.txt', skiprows = 2).T
lamda = R_data[0]
# lnear interpoation prevents interpoation of TRA values outside 0-100%
spectrum_function_list.append( extrap(lamda, R_data[1]/100.0, kind = 'linear' ) )
spectrum_name_list.append('0 deg Reflection')

spectrum_function_list.append( extrap(lamda, R_data[2]/100.0, kind = 'linear' ) )
spectrum_name_list.append('30 deg S-polarization Reflection')
spectrum_function_list.append( extrap(lamda, R_data[3]/100.0, kind = 'linear' ) )
spectrum_name_list.append('30 deg P-polarization Reflection')

spectrum_function_list.append( extrap(lamda, R_data[4]/100.0, kind = 'linear' ) )
spectrum_name_list.append('40 deg S-polarization Reflection')
spectrum_function_list.append( extrap(lamda, R_data[5]/100.0, kind = 'linear' ) )
spectrum_name_list.append('40 deg P-polarization Reflection')

spectrum_function_list.append( extrap(lamda, R_data[6]/100.0, kind = 'linear' ) )
spectrum_name_list.append('50 deg S-polarization Reflection')
spectrum_function_list.append( extrap(lamda, R_data[7]/100.0, kind = 'linear' ) )
spectrum_name_list.append('50 deg P-polarization Reflection')

spectrum_function_list.append( extrap(lamda, R_data[8]/100.0, kind = 'linear' ) )
spectrum_name_list.append('60 deg S-polarization Reflection')
spectrum_function_list.append( extrap(lamda, R_data[9]/100.0, kind = 'linear' ) )
spectrum_name_list.append('60 deg P-polarization Reflection')

spectrum_function_list.append( extrap(lamda, R_data[10]/100.0, kind = 'linear' ) )
spectrum_name_list.append('70 deg S-polarization Reflection')
spectrum_function_list.append( extrap(lamda, R_data[11]/100.0, kind = 'linear' ) )
spectrum_name_list.append('70 deg P-polarization Reflection')


T_data = loadtxt('Transmission_10nm_CuAg_on_silica_substrate.txt', skiprows = 1).T
spectrum_function_list.append( extrap(T_data[0], T_data[1]/100.0, kind = 'linear' ) )
spectrum_name_list.append('0 deg Transmission')


### this function is what we are creating as
def spectrum_list_generator(lamda): # for measured or simulated spectra

	spectrum_list = [spectrum_function(lamda) for spectrum_function in spectrum_function_list] # why not just directly use the spectrum_function_list?
	#because this allows us in the future to use spectra with different Wavelength ranges!
	return spectrum_list #must be a list of spectra  matching the params




#################################### illumination conditions
def generator_function(self,lamda): # layer geometry and neighbor properties

	thickness_list[layer_index_of_fit] = self.thickness
	# order must match the spectrum_list_generator
	param_list = []
	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 0 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 30 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 30 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 40 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 40 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 50 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 50 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )

	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 60 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 60 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 70 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'R'} )


	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 70 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 1.0,
			'spectrum' : 'R'} )



	param_list.append( {
			'lamda' : lamda,
			'snell_angle_front' : 0 * pi/180,
			'layer_index_of_fit' : layer_index_of_fit,
			'nk_f_list' : nk_f_list,
			'thickness_list' : thickness_list,
			'coherency_list' : coherency_list,
			'tm_polarization_fraction' : 0.0,
			'spectrum' : 'T'} )



	return param_list

###### okay so this goofy stuff should allow us to update the thickness of our layer from other places
class updatable_parameter_list_generator:
	def __init__(self):
		self.thickness = 0.0
	def __call__(self,lamda):
		return self.generator_function(self,lamda)


parameter_list_generator = updatable_parameter_list_generator()
parameter_list_generator.thickness = thickness_list[layer_index_of_fit]
parameter_list_generator.generator_function = generator_function





#############################
