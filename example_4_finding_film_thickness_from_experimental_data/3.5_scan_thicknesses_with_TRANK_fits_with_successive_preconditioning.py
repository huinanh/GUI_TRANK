

from TRANK import  try_mkdir,  error_adaptive_iterative_fit_spectra
from numpy import  inf, arange, loadtxt, pi, logspace, log10
from basic_setup import  spectrum_list_generator,   parameter_list_generator, lamda_min, lamda_max




min_thickness = 5
max_thickness = 40

#####uniform spacing
film_thickness_list = arange(min_thickness, max_thickness+.0001 , 1)

####log spacing
#number_of_points = 10
#film_thickness_list = logspace(log10(min_thickness),log10(max_thickness), number_of_points)

### single point thickness
#film_thickness_list = [9.0, 10.5, 11.0, 11.5]


dlamda_min = 5
dlamda_max = 100
delta_weight = 0.1/dlamda_min
lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)

show_plots = False
max_passes = 20 # limit the passes for a quicker scan, 0 means it will guess the amount needed.

if __name__=='__main__':
	if show_plots:
		from matplotlib.pylab import show

	def fit_nk_f(lamda):
		return (1.0 + 0.1j) + 0.0*lamda

	######## KK preconditioning is almost always a good idea, unless you have a metal, which it fails badly for
	film_thickness = film_thickness_list[0]
	print ('\n>>>>>>>> Film Thickness:',film_thickness ,'<<<<<<<<\n')
	parameter_list_generator.thickness = film_thickness
	data_directory = 'TRANK_nk_fit_%f_nm/'%film_thickness
	try_mkdir(data_directory)

	fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
							nk_f_guess = fit_nk_f,
							spectrum_list_generator = spectrum_list_generator,
							parameter_list_generator = parameter_list_generator,
							lamda_min = lamda_min,
							lamda_max = lamda_max,
							dlamda_min = dlamda_min,
							dlamda_max = dlamda_max,
							max_passes = 1,
							delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'linear',
							adaptation_threshold_max = 0.01, adaptation_threshold_min = 0.001,
							use_reducible_error = True,
							method='least_squares',
							KK_compliant = False,
							reuse_mode = False,
							zero_weight_extra_pass = False,
							data_directory = data_directory,
							verbose = True, make_plots = True, show_plots = show_plots,
							nk_spectrum_file_format = 'TRANK_nk_precondition_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_precondition_pass_%i.pdf' )


	for film_thickness in film_thickness_list:
		print ('>>>>>>>> Film Thickness:',film_thickness ,'<<<<<<<<')
		parameter_list_generator.thickness = film_thickness
		data_directory = 'TRANK_nk_fit_%f_nm/'%film_thickness
		try_mkdir(data_directory)

		fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
							nk_f_guess = fit_nk_f,
							spectrum_list_generator = spectrum_list_generator,
							parameter_list_generator = parameter_list_generator,
							lamda_min = lamda_min,
							lamda_max = lamda_max,
							dlamda_min = dlamda_min,
							dlamda_max = dlamda_max,
							max_passes = max_passes,
							lamda_list = lamda_list,## reuse the old lamda points
							delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'cubic',
							adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.0005,
							use_reducible_error = True,
							method='least_squares',
							KK_compliant = False,
							reuse_mode = True,
							zero_weight_extra_pass = False,
							data_directory = data_directory,
							verbose = True, make_plots = True, show_plots = show_plots,
							nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' )
