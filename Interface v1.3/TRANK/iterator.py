from astropy.wcs.docstrings import naxis





def error_adaptive_iterative_fit_spectra(
			nk_f_guess,
			spectrum_list_generator,
			parameter_list_generator,
			lamda_min,
			lamda_max,
			dlamda_min,
			dlamda_max,
			delta_weight = 0.1, tolerance = 1e-5, k_weight_fraction = 1.0,
			adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.0005, adaptation_percentile = 85,
			max_passes = 0,
			lamda_list = [],
			use_reducible_error = True,
			reuse_mode = False,
			KK_compliant = False,
			interpolation_type = 'cubic',
			no_negative = True,
			interpolate_to_fine_grid_at_end = True,
			threads = 0,
			delete_low_error_points = False,
			max_points = 30,
			zero_weight_extra_pass = False, data_directory ='TRANK_nk_fit/', method = 'least_squares', verbose = True, write_nk_files = True,
			make_plots = True, show_plots = True, #nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', 
			#rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' ,
			Gui_mode=False,
			fig=None,
			fig2=None,
			fig3=None,
			TB1=None,
			input_data=None,
			test_setup=None,
			QCoreApplication=None,
			if_e=False
			):


	from TRANK import (fit_spectra_nk_sqr, fit_spectra_nk_sqr_KK_compliant,
						rms_error_spectrum, reducible_rms_error_spectrum, nk_plot, error_plot, try_mkdir)
	from TRANK import compute_coarse_and_fine_grid
	from time import time
	from numpy import floor, log2, ceil, linspace, diff, sqrt, mean, array, savetxt, percentile, argsort
	from copy import deepcopy
	if write_nk_files or make_plots: try_mkdir(data_directory)
    
    
    
    
	
	#if show_plots:
	#	from matplotlib.pylab import show
	error_axes=[] #store axes and present them when called
	nk_axes=[]
	error_ax=None
	nax=None
	kax=None
	fna_rms_spectrum=[] 
	fna_irr_rms_spectrum=[]
	if reuse_mode == False: #picks lambda points accordingly
		lamda_list, lamda_fine = compute_coarse_and_fine_grid(dlamda_max, dlamda_min, lamda_max, lamda_min)
		dlamda_min = lamda_fine[1]-lamda_fine[0]
		dlamda_max = lamda_list[1]-lamda_list[0]

		power_of_2 = int(round( log2(dlamda_max/dlamda_min) ))
		if max_passes == 0:
			passes = int(power_of_2) + 1
		else:
			passes = int(max_passes)



	else:
		lamda_coarse, lamda_fine = compute_coarse_and_fine_grid(dlamda_max, dlamda_min, lamda_max, lamda_min)
		# Determines the fine grid based on the smallest delta lambda you have
		dlamda_min_found = min(diff(lamda_list))
		power_of_2 = int(round( log2(dlamda_min_found/dlamda_min) ))
		#print( log2(dlamda_min_found/dlamda_min)  )
		if False: # in the past it made sense to retrive the dlamda_min and max from data
			dlamda_min = dlamda_min_found/(2**power_of_2) # finds power of two spacing to match your dlamda_min
			dlamda_max = dlamda_max_found = max(diff(lamda_list))

		if max_passes == 0:
			# this part guesses how many passes are required to reach the finest grid level
			passes = max( int(power_of_2) + 1, 2) # this makes sure that it runs on restart!
		else:
			passes = int(max_passes)


	if zero_weight_extra_pass: # this will fail if the num new points conidtion is met
		passes+=1


	fit_nk_f = nk_f_guess
    
    #Print in TextEdit instead of Consol
	#print ('dlamda_max:',dlamda_max )
	TB1.append('dlamda_max:'+str(dlamda_max) )
	#print ('dlamda_min:',dlamda_min )
	TB1.append('dlamda_min:'+str(dlamda_min) )

	# literally jury rigging the conidtion so it starts the loop, ugly, but cleaner than the alternatives
	num_new_points = len(lamda_list)
	total_iteration_time = 0.0
	pass_number = 1
	new_lamda_list = [] # we add no new lamda points for the first pass
	indicies_to_delete = []
	'''
	from new_basic_setup import dyn_basic_setup
	test_setup=dyn_basic_setup(thickness=input_data[0],R_dir=input_data[2],
								T_dir=input_data[3])
	'''
	while pass_number <= passes and num_new_points > 0:

		## delete pointless low error points
		if delete_low_error_points and len(indicies_to_delete)>0:
			#print('Deleted Points:', array(lamda_list)[indicies_to_delete])
			TB1.append('Deleted Points:'+str(array(lamda_list)[indicies_to_delete]))
			for index in indicies_to_delete:
				lamda_list.pop(index)


		## add new lamda points from last pass, does nothing if it is the first pass
		lamda_list = sorted(new_lamda_list+list(lamda_list))
		if pass_number > 1:
			#print('New Points:', new_lamda_list)
			TB1.append('New Points:'+str(new_lamda_list))
			#print('--> Points Added: ', num_new_points)
			TB1.append('--> Points Added: '+str(num_new_points))



		#print('-----------> Pass %i/%i' % (pass_number,passes))
		TB1.append('-----------> Pass %i/%i' % (pass_number,passes))
		#print('--> Fitting %i Points' % len(lamda_list))
		TB1.append('--> Fitting %i Points' % len(lamda_list))

		
		# here we build the inputs for the fitter
		inputs = dict(lamda_list = lamda_list,
					spectrum_list_generator = spectrum_list_generator,
					parameter_list_generator = parameter_list_generator,
					nk_f_guess = fit_nk_f,
					delta_weight = delta_weight,
					tolerance = tolerance,
					no_negative = no_negative,
					k_weight_fraction = k_weight_fraction,
					interpolation_type = interpolation_type, method = method, threads = threads
					,input_data=input_data,
					test_setup=test_setup,
					TB1=TB1,
					)

		t0 = time()
		if KK_compliant:
			inputs.update(dict(lamda_fine = lamda_fine))
			fit_nk_f = fit_spectra_nk_sqr_KK_compliant(**inputs ) # <-----
		else:
			fit_nk_f = fit_spectra_nk_sqr(**inputs)  # <-----
		pass_time = time()-t0

		total_iteration_time += pass_time
		#print('Pass Time: %.1f seconds'%pass_time)
		TB1.append('Pass Time: %.1f seconds'%pass_time)
		#QCoreApplication.processEvents()


		rms_spectrum = rms_error_spectrum(lamda_list = lamda_list,
							nk_f = fit_nk_f,
							spectrum_list_generator = spectrum_list_generator,
							parameter_list_generator = parameter_list_generator, threads = threads,
							input_data=input_data, test_setup=test_setup, get_data=True)
		error_R=[]
		error_T=[]
		temp_rms=[]
		for tu in rms_spectrum:
			temp_rms.append(tu[0])
			error_R.append(tu[1])
			error_T.append(tu[2])
		rms_spectrum=temp_rms
		net_rms = sqrt( mean( array(rms_spectrum)**2 ) )
		max_rms = max(rms_spectrum)

		rms_spectrum_fine = rms_error_spectrum(lamda_list = lamda_fine,
						nk_f = fit_nk_f,
						spectrum_list_generator = spectrum_list_generator,
						parameter_list_generator = parameter_list_generator, threads = threads,
						input_data=input_data, test_setup=test_setup)
		net_rms_fine = sqrt( mean( array(rms_spectrum_fine)**2 ) )


		nk = fit_nk_f(lamda_list)
		if use_reducible_error == False:
			reducible_error_spectrum = []
			adaptation_threshold = max( min(percentile(rms_spectrum,adaptation_percentile),adaptation_threshold_max) , adaptation_threshold_min)
			if write_nk_files: 
				savetxt(data_directory+'fit_nk.txt',array([lamda_list, nk.real, nk.imag, array(rms_spectrum)*100.0]).T)
				fna_rms_spectrum+=(array(rms_spectrum)*100).tolist()
		else:
			reducible_error_spectrum, irreducible_error_spectrum = reducible_rms_error_spectrum(
		 						lamda_list = lamda_list,
								nk_f = fit_nk_f,
								spectrum_list_generator = spectrum_list_generator,
								parameter_list_generator = parameter_list_generator, threads = threads,
								input_data=input_data, test_setup=test_setup)
			adaptation_threshold = max( min(percentile(reducible_error_spectrum, adaptation_percentile),adaptation_threshold_max) , adaptation_threshold_min)
			if write_nk_files: 
				savetxt(data_directory+'fit_nk.txt',array([lamda_list, nk.real, nk.imag, array(rms_spectrum)*100.0, array(reducible_error_spectrum)*100]).T)
				fna_rms_spectrum+=(array(rms_spectrum)*100).tolist()
				fna_irr_rms_spectrum+=(array(irreducible_error_spectrum)*100).tolist()


		#print('Fine Grid Net RMS Error: %f %%' % (net_rms_fine*100))
		TB1.append('Fine Grid Net RMS Error: %f %%' % (net_rms_fine*100))
		#print('--> Net RMS Error: %f %%' % (net_rms*100))
		TB1.append('--> Net RMS Error: %f %%' % (net_rms*100))
		#print('--> Adaptation Threshold: %f %%' % (adaptation_threshold* 100))
		TB1.append('--> Adaptation Threshold: %f %%' % (adaptation_threshold* 100))

		
        
		if make_plots:
			if error_ax:
				fig2.delaxes(error_ax)
			elif fig2.gca()!=None:
				for ax in fig2.axes:
					fig2.delaxes(ax)
			error_ax= error_plot(lamda_list = lamda_list, rms_spectrum = rms_spectrum,fig2=fig2,
							adaptation_threshold = adaptation_threshold,
							adaptation_threshold_min = adaptation_threshold_min,
							adaptation_threshold_max = adaptation_threshold_max,
							reducible_error_spectrum = reducible_error_spectrum,
							#file_name = data_directory + rms_spectrum_file_format % pass_number,
							title_string = 'Thickness= %.3fnm\nPass %i\nNon-Uniform RMS Error = %.3f %%\nUniform Fine RMS Error = %.3f %%' %
							( float(input_data[0]),pass_number, net_rms*100, net_rms_fine*100),
							show_plots = show_plots )
			error_axes.append(error_ax)
			a=len(error_axes)
            
			if nax:
				fig.delaxes(nax)
				fig.delaxes(kax)
			elif fig.gca()!=None:
				for ax in fig.axes:
					fig.delaxes(ax)
			    
            	
			nax, kax= nk_plot(lamda_list = lamda_list, lamda_fine = lamda_fine,fig=fig, nkf = fit_nk_f,
				#file_name = data_directory + nk_spectrum_file_format % pass_number ,
				title_string='Thickness= %.3fnm\nTRANK Pass %i' % (float(input_data[0]),pass_number), 
				show_nodes = True, show_plots = show_plots,if_e=if_e)
			nk_axes.append(nax)
			nk_axes.append(kax)
			
			sim_R=[]
			real_R=[]
			sim_T=[]
			real_T=[]
			for lamda,err_r,err_t in zip(lamda_list,error_R,error_T):
				spectrum=test_setup.spectrum_list_generator(lamda)
				sim_R.append(spectrum[0])
				real_R.append(spectrum[0]+err_r)
				sim_T.append(spectrum[-1])
				real_T.append(spectrum[-1]+err_t)
			if fig3.gca()!=None:
				for ax in fig3.axes:
					fig3.delaxes(ax)
			sim_ax=fig3.add_subplot(121)
			sim_ax_2=fig3.add_subplot(122)
			fig3.tight_layout(pad = 0.8)
			
		
		
			sim_ax.plot(lamda_list, sim_R, linewidth = 2.0, color = 'k', linestyle = '', marker = 'o', markersize = .4, zorder = 1)
			sim_ax.plot(lamda_list, real_R, linewidth = 2.0, color = 'r', label ='n', linestyle = '-', zorder = 0.9)
			
			fig3.subplots_adjust(top = 0.85,
					bottom = 0.25,
					left = 0.10,
					right = 0.95,
					wspace=0.5,
					hspace=0.2)
			sim_ax.set_xlabel('Wavelength (nm)')
			sim_ax.set_ylabel('Reflectance')
			sim_ax.minorticks_on()
			
			sim_ax_2.plot(lamda_list, sim_T, linewidth = 2.0, color = 'k', linestyle = '', marker = 'o', markersize = .4, zorder = 1)
			sim_ax_2.plot(lamda_list, real_T, linewidth = 2.0, color = 'r', label ='n', linestyle = '-', zorder = 0.9)
			
			sim_ax_2.set_xlabel('Wavelength (nm)')
			sim_ax_2.set_ylabel('Transmittance')
			sim_ax_2.minorticks_on()
			fig3.suptitle('Normal Incidence', fontsize=12,x=0.53,y=1.00)
		    
			if show_plots:
				fig.canvas.draw()
				fig2.canvas.draw()
				fig3.canvas.draw()
				#QCoreApplication.processEvents()

		############ adaptation
		if use_reducible_error:
			adaptation_spectrum = reducible_error_spectrum
		else:
			adaptation_spectrum = rms_spectrum


		refinement_method = 'near_worst'
		#### with our adaptation selection method set, we find new points
		if refinement_method == 'near_worst':
			new_lamda_list = []
			for i in range(len(lamda_list)-1):
				if (adaptation_spectrum[i] > adaptation_threshold) or (adaptation_spectrum[i+1] > adaptation_threshold): # should we refine?
					if (lamda_list[i+1] - lamda_list[i]) > dlamda_min: # if the gap is bigger than the minimum, then it is allowed to refine
						new_lamda = (lamda_list[i]+lamda_list[i+1])/2.0
						new_lamda_list.append( new_lamda)

		elif refinement_method == 'interpolate_and_check_all':
			test_lamda_list = []
			for i in range(len(lamda_list)-1):
				if (lamda_list[i+1] - lamda_list[i]) > dlamda_min: # if the gap is bigger than the minimum, then it is allowed to refine
					test_lamda_list.append(  (lamda_list[i]+lamda_list[i+1])/2.0 )

			if use_reducible_error :
				test_error_spectrum, test_irreducible_error_spectrum = reducible_rms_error_spectrum(
			 						lamda_list = test_lamda_list,
									nk_f = fit_nk_f,
									spectrum_list_generator = spectrum_list_generator,
									parameter_list_generator = parameter_list_generator, threads = threads,
									input_data=input_data, test_setup=test_setup)
			else:
				test_error_spectrum = rms_error_spectrum(lamda_list = test_lamda_list,
									nk_f = fit_nk_f,
									spectrum_list_generator = spectrum_list_generator,
									parameter_list_generator = parameter_list_generator, threads = threads,
									input_data=input_data, test_setup=test_setup)
			#sorted_indices =  argsort(test_error_spectrum)
			new_lamda_list = []
			for i in range(len(test_lamda_list)):
				if (test_error_spectrum[i] > adaptation_threshold) :
					new_lamda_list.append( test_lamda_list[i])

		#### we combine the new points with the old at the start of the next pass
		### this is important to have here for the termination condition
		num_new_points = len(new_lamda_list)

		############ adaptation
		if delete_low_error_points:
			if ( (num_new_points + len(lamda_list)) > max_points):
				n_delete = num_new_points+len(lamda_list) - max_points
				sorted_indices =  argsort(adaptation_spectrum)
				### remove edge indices
				sorted_indices_without_edge_values = list(sorted_indices)
				sorted_indices_without_edge_values.remove(0)
				sorted_indices_without_edge_values.remove(len(adaptation_spectrum)-1)

				# now we remove any that would make a gap that is too large
				for index_index in range(len(sorted_indices_without_edge_values)-1,-1,-1):
					index_to_check = sorted_indices_without_edge_values[index_index]
					if (lamda_list[index_to_check+1] - lamda_list[index_to_check-1]) > dlamda_max:
						 del sorted_indices_without_edge_values[index_index] # we can't consider it, would make a large gap

				indicies_to_delete  = sorted_indices_without_edge_values[0:n_delete]
				indicies_to_delete.sort(reverse = True)



		#### doing the stuff for the last extra pass if there is one
		if zero_weight_extra_pass:
			if (pass_number +1) == passes: # normal zero_weight_extra_pass , just finished second to last pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number += 1

			elif num_new_points == 0 and pass_number < passes: # test if terminates early, but still needs that extra pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number = passes # skip to last passes
				#print('--> Skipping to extra pass due to early conidtion statisfaction')
				TB1.append('--> Skipping to extra pass due to early conidtion statisfaction')
		else:
			pass_number += 1

	TB1.append('Total Iterating Time: %.1f seconds'%total_iteration_time)

	

    
	if interpolate_to_fine_grid_at_end and write_nk_files:
		TB1.append('Interpolating to fine grid and saving...')
		#QCoreApplication.processEvents()
		nk = fit_nk_f(lamda_fine)

		if use_reducible_error == False:
			savetxt(data_directory+'fit_nk_fine.txt',array([lamda_fine, nk.real, nk.imag, array(rms_spectrum_fine)*100.0]).T)
		else:
			reducible_error_spectrum_fine, irreducible_error_spectrum = reducible_rms_error_spectrum(
		 						lamda_list = lamda_fine,
								nk_f = fit_nk_f,
								spectrum_list_generator = spectrum_list_generator,
								parameter_list_generator = parameter_list_generator, threads = threads,
								input_data=input_data, test_setup=test_setup)
			savetxt(data_directory+'fit_nk_fine.txt',array([lamda_fine, nk.real, nk.imag, array(rms_spectrum_fine)*100.0,  array(reducible_error_spectrum_fine)*100.0 ]).T)

	TB1.append('Completed!')
	
	if Gui_mode==1:
		return error_axes, nk_axes
	elif Gui_mode==2:
		return fna_rms_spectrum, fit_nk_f, lamda_list
	elif Gui_mode==3:
		return fna_irr_rms_spectrum, fit_nk_f, lamda_list
	elif Gui_mode==4:
		return error_axes, nk_axes,fit_nk_f, lamda_list
	
	else:
		return fit_nk_f, lamda_list














def error_adaptive_iterative_fit(
			nk_f_guess,
			TR_pair_list_generator,
			parameter_list_generator,
			lamda_min,
			lamda_max,
			dlamda_min,
			dlamda_max,
			delta_weight = 0.1, tolerance = 1e-5,
			adaptation_threshold_max = 0.01, adaptation_threshold_min = 0.0005,
			max_passes = 0,
			extra_passes = 0,
			lamda_list = [],
			use_reducible_error = True,
			reuse_mode = False,
			KK_compliant = False,
			interpolation_type = 'cubic',
			zero_weight_extra_pass = False, data_directory ='TRANK_nk_fit/', method = 'least_squares', verbose = True,
			make_plots = True, show_plots = True, nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf'  ):


	from TRANK import (fit_TRA_nk_sqr, fit_TRA_nk_sqr_KK_compliant,
						rms_TRA_error_spectrum, reducible_rms_TRA_error_spectrum, nk_plot, error_plot, try_mkdir)
	from time import time
	from numpy import floor, log2, ceil, linspace, diff, sqrt, mean, array, savetxt, percentile
	try_mkdir(data_directory)

	if show_plots:
		from matplotlib.pylab import show

	if reuse_mode == False: #picks lambda points accordingly
		ncoarse = ceil((lamda_max - lamda_min)/dlamda_max)
		dlamda_max =  (lamda_max - lamda_min)/ncoarse
		lamda_list = linspace(lamda_min,  lamda_max, ncoarse+1)
		#print(lamda_list)
		power_of_2 = int(round( log2(dlamda_max/dlamda_min) ))
		#print(log2(dlamda_max/dlamda_min), power_of_2)
		dlamda_min = dlamda_max/(2**power_of_2)
		lamda_fine = linspace(lamda_min,  lamda_max, ncoarse*(2**power_of_2)+1)
		if max_passes == 0:
			passes = int(power_of_2) + 1
		else:
			passes = int(max_passes)

	else:
		dlamda_min_found = min(diff(lamda_list))
		power_of_2 = int(round( log2(dlamda_min_found/dlamda_min) ))
		#print( log2(dlamda_min_found/dlamda_min)  )
		dlamda_min = dlamda_min_found/(2**power_of_2)
		#print ('dlamda_min', dlamda_min)
		nfine = ceil((lamda_max - lamda_min)/dlamda_min)
		#print ('nfine', nfine)
		lamda_fine = linspace(lamda_min,  lamda_max, nfine+1)
		#print ('lamda_fine', lamda_fine)

		if max_passes == 0:
			passes = max( int(power_of_2) + 1, 2) # this makes sure that it runs on restart!
		else:
			passes = int(max_passes)



	if zero_weight_extra_pass: # this will fail if the num new points conidtion is met
		passes+=1


	fit_nk_f = nk_f_guess

	print ('dlamda_max:',dlamda_max )
	print ('dlamda_min:',dlamda_min )

	num_new_points = len(lamda_list)
	total_iteration_time = 0.0
	pass_number = 1
	while pass_number <= passes and num_new_points > 0:

		print('-----------> Pass %i/%i' % (pass_number,passes))
		print('--> Fitting %i Points' % len(lamda_list))


		# here we build the inputs for the fitter
		inputs = dict(lamda_list = lamda_list,
					TR_pair_list_generator = TR_pair_list_generator,
					parameter_list_generator = parameter_list_generator,
					nk_f_guess = fit_nk_f,
					delta_weight = delta_weight,
					tolerance = tolerance,
					interpolation_type = interpolation_type, method = method)

		t0 = time()
		if KK_compliant:
			inputs.update(dict(lamda_fine = lamda_fine))
			fit_nk_f = fit_TRA_nk_sqr_KK_compliant(**inputs ) # <-----
		else:
			fit_nk_f = fit_TRA_nk_sqr(**inputs)  # <-----
		pass_time = time()-t0

		total_iteration_time += pass_time
		print('Pass Time: %.1f seconds'%pass_time)


		rms_spectrum = rms_TRA_error_spectrum(lamda_list = lamda_list,
							nk_f = fit_nk_f,
							TR_pair_list_generator = TR_pair_list_generator,
							parameter_list_generator = parameter_list_generator)
		net_rms = sqrt( mean( array(rms_spectrum)**2 ) )
		max_rms = max(rms_spectrum)

		rms_spectrum_fine = rms_TRA_error_spectrum(lamda_list = lamda_fine,
						nk_f = fit_nk_f,
						TR_pair_list_generator = TR_pair_list_generator,
						parameter_list_generator = parameter_list_generator)
		net_rms_fine = sqrt( mean( array(rms_spectrum_fine)**2 ) )

		### saving the pass data
		nk = fit_nk_f(lamda_list)
		savetxt(data_directory+'fit_nk.txt',array([lamda_list, nk.real, nk.imag, array(rms_spectrum)*100.0]).T)

		if use_reducible_error:
			reducible_error_spectrum, irreducible_error_spectrum = reducible_rms_TRA_error_spectrum(
		 						lamda_list = lamda_list,
								nk_f = fit_nk_f,
								TR_pair_list_generator = TR_pair_list_generator,
								parameter_list_generator = parameter_list_generator)
			adaptation_threshold = max( min(percentile(reducible_error_spectrum,85),adaptation_threshold_max) , adaptation_threshold_min)
		else:
			reducible_error_spectrum = []
			adaptation_threshold = max( min(percentile(rms_spectrum,85),adaptation_threshold_max) , adaptation_threshold_min)


		print('Fine Grid Net RMS Error: %f %%' % (net_rms_fine*100))
		print('--> Net RMS Error: %f %%' % (net_rms*100))
		print('--> Adaptation Threshold: %f %%' % (adaptation_threshold* 100))


		if make_plots:
			err_fig = error_plot(lamda_list = lamda_list, rms_spectrum = rms_spectrum,
							adaptation_threshold = adaptation_threshold,
							adaptation_threshold_min = adaptation_threshold_min,
							adaptation_threshold_max = adaptation_threshold_max,
							reducible_error_spectrum = reducible_error_spectrum,
							file_name = data_directory+ rms_spectrum_file_format % pass_number,
							title_string = 'Pass %i: Net RMS Error = %.3f %%' %( pass_number, net_rms*100),
							show_plots = show_plots )

			nk_fig = nk_plot(lamda_list = lamda_list, lamda_fine = lamda_fine, nkf = fit_nk_f,
				file_name = data_directory + nk_spectrum_file_format % pass_number ,title_string='TRANK Pass %i' % pass_number, show_nodes = True, show_plots = show_plots)

			if show_plots:
				show()


		if use_reducible_error:
			adaptation_spectrum = reducible_error_spectrum
		else:
			adaptation_spectrum = rms_spectrum
		############ adaptation
		new_lamda_list = []
		#adaptation_threshold = max(rms_spectrum )/2.0
		for i in range(len(lamda_list)-1):
			if (adaptation_spectrum[i] > adaptation_threshold) or (adaptation_spectrum[i+1] > adaptation_threshold): # should we refine?
				if (lamda_list[i+1] - lamda_list[i]) > dlamda_min: # if the gap is bigger than the minimum, then it is allowed to refine
					new_lamda = (lamda_list[i]+lamda_list[i+1])/2.0
					new_lamda_list.append( new_lamda)

		#### now we combine the new points with the old
		num_new_points = len(new_lamda_list)
		print('New Points:', new_lamda_list)
		print('--> Points Added: ', num_new_points)

		lamda_list = sorted(new_lamda_list+list(lamda_list))

		#### doing the stuff for the last extra pass if there is one
		if zero_weight_extra_pass:
			if (pass_number +1) == passes: # normal zero_weight_extra_pass , just finished second to last pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number += 1

			elif num_new_points == 0 and pass_number < passes: # test if terminates early, but still needs that extra pass
				delta_weight = 0.0
				tolerance = 1e-8
				num_new_points = 1 # jury rig it so it continues regardless of state of convergence
				pass_number = passes # skip to last passes
				print('--> Skipping to zero weight extra pass due to early conidtion statisfaction')
		else:
			pass_number += 1




	print('Total Iterating Time: %.1f seconds'%total_iteration_time)
	nk = fit_nk_f(lamda_fine)
	savetxt(data_directory+'fit_nk_fine.txt',array([lamda_fine, nk.real, nk.imag, array(rms_spectrum_fine)*100.0]).T)
    	
	return fit_nk_f
