###########################################################




show_plots = True



########### get stucture
from basic_setup  import  spectrum_list_generator,   parameter_list_generator, spectrum_name_list


data_directory = 'TRANK_nk_fit/'
from TRANK import functionize_nk_file, TMM_spectra, error_plot, reducible_rms_error_spectrum, rms_error_spectrum


if __name__=='__main__':


	from numpy import loadtxt, array
	lamda_fine = loadtxt(data_directory+'fit_nk_fine.txt', usecols = [0], unpack = True)
	lamda_list = loadtxt(data_directory+'fit_nk.txt', usecols = [0], unpack = True)
	fit_nk_f = functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'cubic')



	spectra_from_fit = array(TMM_spectra(lamda_list = lamda_fine,
										nk_f = fit_nk_f,
										parameter_list_generator = parameter_list_generator))
	#spectra_from_fit[lamda][ param spectrum] # computational reasons why this is ordered is this way


	original_spectra = spectrum_list_generator(lamda_fine)
	#original_spectra[spectrum number][lamda] # again, computational reasons why this is ordered is this way


	from matplotlib.pylab import figure, plot, legend, minorticks_on, gca, title, tight_layout, ylabel, xlabel, savefig, show

	# you could do this if you had a mixed bag of conifurations
	#num_parameters = len(parameter_list_generator(lamda = lamda_fine[0]))
	#for parameter_number in range(num_parameters):
	figure(figsize =(3.2 * 2, 2.5 * 2), dpi = 220*2.0/3.0 )

	for spectrum_index in range(len(spectrum_name_list)):


		name = spectrum_name_list[spectrum_index]

		fit_data = spectra_from_fit[:, spectrum_index]
		exp_data = original_spectra[spectrum_index]

		p = plot(lamda_fine, fit_data*100.0, label = 'Fit '+name,            linewidth = 1.0, linestyle = '--')

		plot(   lamda_fine, exp_data *100.0 , label = 'Exp. '+name, color = p[0].get_color(), linewidth = 1.0)


	legend(loc= 'center right', labelspacing = 0.09, handletextpad = 0.3, handlelength = 0.8, fontsize = 8 )
	minorticks_on()
	gca().tick_params(axis='both', which='major', labelsize=10)

	ylabel('(%)',fontsize =10)
	xlabel('Wavelength (nm)',fontsize =10)
	#angle = angle_list[angle_index]
	#title('%i$^{\circ}$'%angle, fontsize =12)

	gca().set_ylim(0.0,100.0)
	gca().set_xlim(min(lamda_fine),max(lamda_fine))

	tight_layout(pad=0.1)
	savefig(data_directory+'spectra_comparison.pdf' ,dpi = 600, transparent = True)

	########
	rms_spectrum = rms_error_spectrum(lamda_list = lamda_list, nk_f = fit_nk_f,
					spectrum_list_generator = spectrum_list_generator,
					parameter_list_generator = parameter_list_generator)

	reducible_error_spectrum, irreducible_error_spectrum = reducible_rms_error_spectrum(lamda_list= lamda_list, nk_f = fit_nk_f,
								spectrum_list_generator = spectrum_list_generator,
								parameter_list_generator = parameter_list_generator)

	rms_spectrum_fine = rms_error_spectrum(lamda_list = lamda_fine, nk_f = fit_nk_f,
						spectrum_list_generator = spectrum_list_generator,
						parameter_list_generator = parameter_list_generator)

	reducible_error_spectrum_fine, irreducible_error_spectrum_fine = reducible_rms_error_spectrum(lamda_list= lamda_fine, nk_f = fit_nk_f,
									spectrum_list_generator = spectrum_list_generator,
									parameter_list_generator = parameter_list_generator)

	error_plot(lamda_list = lamda_list, rms_spectrum = rms_spectrum,
					adaptation_threshold = -1.0, adaptation_threshold_min = -1.0, adaptation_threshold_max = -1.0,
					reducible_error_spectrum = reducible_error_spectrum,
					lamda_fine = lamda_fine, rms_spectrum_fine = rms_spectrum_fine, reducible_error_spectrum_fine = reducible_error_spectrum_fine,
					title_string = '',
					file_name = data_directory+'final_error_spectrum_fine.pdf', zoom_window = [], show_plots = True )

	if show_plots: show()
