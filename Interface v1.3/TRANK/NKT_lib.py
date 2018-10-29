

from tmm import inc_tmm
from numpy import sqrt


def NKT_TMM_spectrum_wrapper(nk_fit, thickness, lamda, snell_angle_front, layer_index_of_fit,  nk_f_list,
							thickness_list, coherency_list, tm_polarization_fraction, spectrum):
	#this is just a fancy wrapper for inc_tmm
	# does the order matter?

	nk_list = [ nk_f(lamda) for nk_f in nk_f_list]
	nk_list[layer_index_of_fit] = nk_fit # overwrite input nk_f with the fit one. basically, it would make the code much uglier if I made an excption

	local_thickness_list = [layer_thickness for layer_thickness in thickness_list]
	local_thickness_list[layer_index_of_fit] = thickness

	te_result = inc_tmm('s', nk_list, local_thickness_list, coherency_list, snell_angle_front, lamda)
	tm_result = inc_tmm('p', nk_list, local_thickness_list, coherency_list, snell_angle_front, lamda)

	T = tm_polarization_fraction * tm_result['T'] + (1.0-tm_polarization_fraction) * te_result['T']
	R = tm_polarization_fraction * tm_result['R'] + (1.0-tm_polarization_fraction) * te_result['R']
	A = 1 - T - R

	if callable(spectrum)==False:
		result_dict = {	'T': T,
						'R': R,
						'A': A }
		result = result_dict[spectrum]
	else:
		result = spectrum(T,R) # allows you to create things like extiction where the spectrum is 1-T

	return result


def NKT_spectrum_TMM_lamda(params): # This has to be at the top level because map is strange and wont pickle onless it is at the top level
	'''Returns a list of the spectrum values for each single point at a Wavelength'''
	lamda = params[0]
	nk    = params[1]
	thickness = params[2]
	parameter_list_generator = params[3]

	list_of_parameters = parameter_list_generator(lamda)
	spectrum_list = []
	for parameters in list_of_parameters :
		spectrum_list.append(  TMM_spectrum_wrapper(nk_fit = nk, thickness = thickness,  **parameters) )
	return spectrum_list


def NKT_spectrum_lamda_error(params): # This has to be at the top level because map is strange and wont pickle onless it is at the top level
	lamda = params[0] # these params are per Wavelength, we could get more ganular parallelism with this!
	nk    = params[1]
	thickness = params[2]
	spectrum_list_generator = params[3]
	parameter_list_generator = params[4]

	spectrum_list = spectrum_list_generator(lamda)
	sqrt_point_multiplicity = sqrt(len(spectrum_list))
	list_of_parameters = parameter_list_generator(lamda)
	sub_error_list=[]
	for spectrum, parameters in zip (spectrum_list, list_of_parameters ):
		spectrum_calculated = NKT_TMM_spectrum_wrapper(nk_fit = nk, thickness = thickness,  **parameters) # these parameters are per measurement, set of parameters for a single point tmm model
		error = (spectrum_calculated - spectrum)/sqrt_point_multiplicity # in the formulation this gets squared later, and we look at the per Wavelength rms error and normalizing it gives portability to weights
		#we'll put wieghting somewhere else so that the formulation is clear
		sub_error_list.append(  error)

	## i should try this later
	#def error( value_and_model_parameters_tuple ):
	#	tmm_parameters = value_and_model_parameters_tuple[1]
	#	value = value_and_model_parameters_tuple[0]
	#	error = (TMM_spectrum(nk_fit = nk,  **tmm_parameters ) - value)/sqrt_point_multiplicity
	#	return  error

	#sub_error_list = map(error, zip(spectrum_list, list_of_parameters ) )

	return sub_error_list


def NKT_TMM_spectra(lamda_list, nk_f,  parameter_list_generator):
	''''''
	#returns data in block like spectra[lamda][ spectrum] there are computational reasons why this order is this way
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	spectra = my_pool.map(NKT_spectrum_TMM_lamda, muh_inputs)
	my_pool.terminate()
	#
	return spectra



def NKT_pointwise_rms_error_sum_wrapper(params): # for each Wavelength, wraps the previous function TRA_lamda_error to quikcly compute error spectrum
		from numpy import  sqrt
		point_error_list = NKT_spectrum_lamda_error(params)
		sum_err_square = 0.0
		for err in point_error_list:
			sum_err_square += err**2
		return sqrt(sum_err_square)



def NKT_rms_error_spectrum(lamda_list, nk_f, thickness, spectrum_list_generator, parameter_list_generator):
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), thickness, spectrum_list_generator, parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	error_spectrum = my_pool.map(NKT_pointwise_rms_error_sum_wrapper, muh_inputs)

	my_pool.terminate()
	return error_spectrum

def NKT_sqr_rms_gradient_at_lamda(params, h_nk = 1e-6):
	from numpy import array
	def shift_nk(params, shift): # i use this code to make it clearer
		return [params[0], params[1]+shift, params[2], params[3], params[4]]

	e_p_0 = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk ))**2
	e_m_0 = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk ))**2

	e_0_p = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*1.0j ))**2
	e_0_m = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk*1.0j ))**2

	e_n = (e_p_0 - e_m_0)/(2.0*h_nk)
	e_k = (e_0_p - e_0_m)/(2.0*h_nk)
	gradient = array([e_n, e_k])
	return gradient


def NKT_sqr_rms_hessian_at_lamda(params, h_nk = 1e-4):
	from numpy import array
	def shift_nk(params, shift): # i use this code to make it clearer
		return [params[0], params[1]+shift, params[2], params[3], params[4]]
	e_0_0 = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  0.0 ))**2

	e_p_0 = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0+0.0j) ))**2
	e_m_0 = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk*( 1.0+0.0j) ))**2

	e_0_p = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 0.0+1.0j) ))**2
	e_0_m = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk*( 0.0+1.0j) ))**2

	e_p_p = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0+1.0j) ))**2
	e_m_m = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0-1.0j) ))**2

	e_p_m = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0-1.0j) ))**2
	e_m_p = NKT_pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0+1.0j) ))**2

	e_nn = (e_p_0 - 2.0*e_0_0 + e_m_0)/(h_nk**2.0)
	e_kk = (e_0_p - 2.0*e_0_0 + e_0_m)/(h_nk**2.0)

	e_nk = (e_p_p - e_p_0 - e_0_p + 2.0*e_0_0 - e_m_0 - e_0_m + e_m_m)/(2.0*h_nk*h_nk) # two difference methods..., the first is supposed to be cheaper
	e_nk = (e_p_p - e_p_m - e_m_p + e_m_m)/(4.0*h_nk*h_nk)

	hessian = array([[e_nn, e_nk],[e_nk, e_kk]])
	return hessian

def NKT_pointwise_reducible_rms_error_sum_wrapper(params): # for each Wavelength, wraps the previous function TRA_lamda_error to quikcly compute error spectrum
	## calculates numeric derivatives
	#does some hessians and matrix match
	#from wikipedia

	e_0_0 = NKT_pointwise_rms_error_sum_wrapper(params = params)**2
	grad = NKT_sqr_rms_gradient_at_lamda(params = params)# h_nk = h_nk)
	hessian = NKT_sqr_rms_hessian_at_lamda(params = params)#, h_nk = h_nk)

	from numpy.linalg import det, inv
	from numpy import dot, array, sqrt

	if det(hessian) > 0.0: # minima predicted
		hessian_inv = inv(hessian)
		#reducible_error = 1/2.0 * dot(grad, dot( hessian_inv, grad))
		#if reducible_error > e_0_0: # can't go negative this way
		#	reducible_error = e_0_0

		#irreducible_error = e_0_0 - reducible_error
		S_change = 1.0/2.0 * dot(grad, dot( hessian_inv, grad))
		Smin = e_0_0 - S_change
		irreducible_error = Smin
		reducible_error = e_0_0 - irreducible_error

		if irreducible_error < 0.0: # removing this might make better results even if unrealisitc
			reducible_error =  e_0_0 # bascially its predicting negative error is posible, nah just meas it can go to zero
			irreducible_error = 0.0

		#if reducible_error < 0.0:
		#	print( params[0], det(hessian), S_change, [reducible_error, irreducible_error])
	else: # what else should i do here, I assume it means the band is caught on a peak?
		irreducible_error = e_0_0
		reducible_error = 0.0
		#print( params[0], det(hessian), [reducible_error, irreducible_error])

	return [reducible_error, irreducible_error]

def NKT_reducible_rms_error_spectrum(lamda_list, nk_f, thickness, spectrum_list_generator, parameter_list_generator):
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), thickness, spectrum_list_generator, parameter_list_generator ) )

	from numpy import array, sqrt
	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	error_spectra = array(my_pool.map(NKT_pointwise_reducible_rms_error_sum_wrapper, muh_inputs)).T
	my_pool.terminate()

	reducible_error_spectrum = sqrt(error_spectra[0])
	irreducible_error_spectrum = sqrt(error_spectra[1])
	return reducible_error_spectrum, irreducible_error_spectrum



def NKT_single_lamda_rms_error_map(lamda, nlist, klist, thickness,  spectrum_list_generator, parameter_list_generator):
	muh_inputs = []
	for n in nlist:
		for k in klist:
			muh_inputs.append( (lamda, n+1.0j*k, thickness, spectrum_list_generator, parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())

	from numpy import reshape, array
	error_list = my_pool.map(NKT_pointwise_rms_error_sum_wrapper, muh_inputs)
	error_map = reshape( array(error_list), (len(nlist), len(klist) ))
	my_pool.terminate()
	return error_map #[nindex,kindex]



###################
def NKT_fit_spectra_nk_sqr(lamda_list, spectrum_list_generator, parameter_list_generator,  nk_f_guess, thickness_guess, delta_weight = 0.1, tolerance = 1e-5, no_negative = True, interpolation_type = 'cubic', method = 'least_squares'):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi, exp, abs, sqrt, array, zeros, savetxt, inf, diff, ones
	from scipy.optimize import root, least_squares, minimize
	from TRANK import extrap_c


	#point_multiplicity = len(spectrum_list_generator(lamda_list[0]))
	#print(point_multiplicity)

	#point_multiplicity_list = [len(spectrum_list_generator(lamda)) for lamda in lamda_list ]
	#point_multiplicity = point_multiplicity_list[0]

	abs_delta_weight = sqrt(delta_weight**2  * (len(lamda_list)/(len(lamda_list)-1.0)))

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	#my_pool = Pool(1)


	def F_error(nk_t_list): # t goes at the end

		c_nk_list = []
		muh_inputs = []
		thickness = nk_t_list[len(lamda_list)]
		for i in range(len(lamda_list)):
			nk = nk_t_list[2*i] + 1.0j*nk_t_list[2*i+1]
			c_nk_list.append(nk)
			muh_inputs.append( (lamda_list[i], nk, thickness, spectrum_list_generator, parameter_list_generator ) )


		error_list_lists = my_pool.map(NKT_spectrum_lamda_error, muh_inputs)

		#combine the sub error lists into
		error_list = []
		for sub_error_list in error_list_lists:
			error_list = error_list + sub_error_list

		delta_array = diff(c_nk_list)*abs_delta_weight

		error_list = error_list + list(delta_array.real) + list(delta_array.imag)

		return error_list

	####### now for a guess list ##############

	nk_t_guess_list = []
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		if no_negative:
			nk_t_guess_list.append(abs(nk.real))
			nk_t_guess_list.append(abs(nk.imag))
		else:
			nk_t_guess_list.append(nk.real)
			nk_t_guess_list.append(nk.imag)
	nk_t_guess_list.append(thickness_guess)

	######### test
	if False:
		print(F_error(nk_t_guess_list))
		print( 0.5*sum(array(F_error(nk_t_guess_list))**2))
	############ nk guessing over, time for creating and minimizing error function
	if method == 'least_squares':
		inputs = dict(fun = F_error,
					x0 = nk_t_guess_list,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2)
		if no_negative:
			inputs.update(dict( bounds = [     [0]*(2*len(lamda_list)) + [0] , [inf] *(2*len(lamda_list)) + [inf]  ]))
		else:
			inputs.update(dict( bounds = [  [-inf]*(2*len(lamda_list)) + [0] , [inf] *(2*len(lamda_list)) + [inf]  ]))
		solution = least_squares(**inputs ).x

	elif method == 'L-BFGS-B':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = nk_t_guess_list,
				method = 'L-BFGS-B',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(   0,inf)] + [(0,inf)] ))
		else:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(-inf,inf)] + [(0,inf)] ))
		solution = minimize(**inputs ).x

	elif method == 'SLSQP':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = nk_t_guess_list,
				method = 'SLSQP',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(   0,inf)] + [(0,inf)] ))
		else:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(-inf,inf)] + [(0,inf)] ))
		solution = minimize(**inputs ).x

	elif method == 'TNC':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = nk_t_guess_list,
				method = 'TNC',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(   0,inf)] + [(0,inf)] ))
		else:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(-inf,inf)] + [(0,inf)] ))
		solution = minimize(**inputs ).x

	else:
		raise ValueError("Invalid minimization method!")

	my_pool.terminate()
	my_pool.close()


	nk_list=[]
	for i in range(len(lamda_list)):
		nk_list.append(solution[2*i] + 1.0j*solution[2*i+1]  )
	thickness = solution[len(lamda_list)]

	fit_nk_f = extrap_c(lamda_list, nk_list, kind = interpolation_type)


	return fit_nk_f, thickness





def KK_lamda(lamda_list, lamda_fine,  k,  cshift = 1e-4) :
	#KK transform in Wavelength
	# cshift is basically scaled so its universal i think
	# testing shows that cubic interpoltion and trapezoind rule excede the accurcy of just using the fine grid due to error cancelation
	from scipy.integrate import simps, trapz
	from scipy.interpolate import griddata

	from numpy import array, zeros, pi
	#print (lamda_list, k)
	#rint (len(lamda_list), len(k))
	k_fine = griddata(array(lamda_list), array(k), (array(lamda_fine)), method='cubic', fill_value = 0.0)
	#print (k_fine)
	n_out = zeros(len(lamda_list)) ### we use the fine lamda grid for the KK integral but we only need to evalute it at the coarse lamba grid points!
	k_over_lamda = k_fine/lamda_fine
	for i in range(len(lamda_list)): #parralelize this in the future!
		lamda_out = lamda_list[i]

		#n_out[i] = 1.0 + 2.0/pi * simps( k_over_lamda * (1.0/ ( (lamda_out/lamda_fine)**2 - 1.0 + cshift*1.0j)), lamda_fine).real
		n_out[i] = 1.0 + 2.0/pi * trapz( k_over_lamda * (1.0/ ( (lamda_out/lamda_fine)**2 - 1.0 + cshift*1.0j)), lamda_fine).real

	return n_out




def NKT_fit_spectra_nk_sqr_KK_compliant(lamda_list, lamda_fine, spectrum_list_generator, parameter_list_generator,  nk_f_guess, thickness_guess,
								delta_weight = 0.1, tolerance = 1e-5, no_negative = True, interpolation_type = 'cubic', method = 'least_squares'):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi,exp,abs,sqrt, array, matmul, loadtxt, zeros, savetxt, inf, diff, ones, mean
	from scipy.optimize import root, least_squares, minimize
	from TRANK import extrap_c


	#point_multiplicity = len(TR_pair_list_generator(lamda_list[0]))
	#print(point_multiplicity)

	abs_delta_weight = sqrt(delta_weight**2  * (len(lamda_list)/(len(lamda_list)-1.0)))

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	#my_pool = Pool(1)


	def F_error(k_p_t_list):
		# the last value is the principle value
		#FYI -> k = array(k_and_p_list[0:-1])
		k = k_p_t_list[0:len(lamda_list)]
		p = k_p_t_list[len(lamda_list)]
		thickness = k_p_t_list[len(lamda_list)+1]
		n = p + KK_lamda(lamda_list = lamda_list, lamda_fine = lamda_fine,  k = k )



		muh_inputs = []
		for i in range(len(lamda_list)):
			nk = n[i]+1.0j*k[i]# double check this works properly later

			muh_inputs.append( (lamda_list[i], nk, thickness, spectrum_list_generator, parameter_list_generator ) )

		#print (zip(lamda_list, c_nk_list))
		error_list_lists = my_pool.map(spectrum_lamda_error, muh_inputs)
		#error_list_lists =my_pool.map(lamda_error, zip(lamda_list, c_nk_list))
		#print (error_list_lists)

		error_list = []
		for sub_error_list in error_list_lists:
			error_list = error_list + sub_error_list


		error_list = error_list + list( abs_delta_weight*diff(n) )   + list( abs_delta_weight * diff(k))

		return error_list

	####### now for a guess list ##############

	guess_k_p_t_list= []
	p = 0.0
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		if no_negative and nk.imag < 0.0:
			k = 0
		else:
			k = nk.imag

		guess_k_p_t_list.append( k)
		p+= nk.real
	# now we put p at the end
	p = p/len(lamda_list) - 1.0   # this is a guess for the principle value
	print ('principle value guess:',p)
	guess_k_p_t_list.append(p)
	guess_k_p_t_list.append(thickness)

	######### test
	if False: print(F_error(guess_k_p_t_list)) #use this to see if the TR_error works
	############ nk guessing over, time for creating and minimizing error function

	if method == 'least_squares':
		inputs = dict(fun = F_error,
					x0 =  guess_k_p_t_list,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2)

		if no_negative:
			inputs.update(dict( bounds = [     [0]*len(lamda_list) + [0,   0] , [inf]*len(lamda_list) + [inf, inf]  ]))
		else:
			inputs.update(dict( bounds = [  [-inf]*len(lamda_list) + [-inf,0] , [inf]*len(lamda_list) + [inf, inf]  ]))
		#if no_negative:
		#	inputs.update(dict( bounds = [len(lamda_list)*[0.0]+[-inf],len(lamda_list)*[inf]+[inf]] ))
		solution = least_squares(**inputs ).x

	elif method == 'SLSQP':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = guess_k_and_p_list,
				method = 'SLSQP',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(   0,inf)] + [(0,inf)] ))
		else:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(-inf,inf)] + [(0,inf)] ))
		solution = minimize(**inputs ).x

	elif method == 'L-BFGS-B':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = guess_k_and_p_list,
				method = 'L-BFGS-B',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = len(lamda_list)*[(0,inf)]  + [(-inf,inf)] ))
		solution = minimize(**inputs ).x

	elif method == 'TNC':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = guess_k_and_p_list,
				method = 'TNC',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = len(lamda_list)*[(0,inf)]  + [(-inf,inf)] ))
		solution = minimize(**inputs ).x

	else:
		raise ValueError("Invalid minimization method!")

	my_pool.terminate()
	my_pool.close()

	k_and_p_list = solution
	k = k_and_p_list[0:-1]
	p = k_and_p_list[-1]
	n = p + KK_lamda(lamda_list = lamda_list, lamda_fine = lamda_fine,  k = k )
	print ('Final principle value:',p)


	fit_nk_f = extrap_c(lamda_list, n + 1.0j*k, kind = interpolation_type)

	return fit_nk_f
