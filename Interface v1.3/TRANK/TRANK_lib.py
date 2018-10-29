# -*- coding: utf-8 -*- 

from tmm import inc_tmm
from numpy import sqrt
from astropy.wcs.docstrings import spec


from _functools import partial





def TMM_spectrum_wrapper(nk_fit, lamda, snell_angle_front, layer_index_of_fit,  nk_f_list,  thickness_list, coherency_list, tm_polarization_fraction, spectrum, get_data=False): #this is just a fancy wrapper for inc_tmm

	nk_list = [ nk_f(lamda) for nk_f in nk_f_list]
	nk_list[layer_index_of_fit] = nk_fit # overwrite input nk_f with the fit one. basically, it would make the code much uglier if I made an excption

	te_result = inc_tmm('s', nk_list, thickness_list, coherency_list, snell_angle_front, lamda)
	tm_result = inc_tmm('p', nk_list, thickness_list, coherency_list, snell_angle_front, lamda)

	T = tm_polarization_fraction * tm_result['T'] + (1.0-tm_polarization_fraction) * te_result['T']
	R = tm_polarization_fraction * tm_result['R'] + (1.0-tm_polarization_fraction) * te_result['R']


	if callable(spectrum)==False:
		A = 1 - T - R
		result_dict = {	'T': T,
						'R': R,
						'A': A }
		result = result_dict[spectrum]
	else:
		result = spectrum(T,R) # allows you to create things like extiction where the spectrum is 1-T

	return result


def spectrum_lamda_error(params): # This has to be at the top level because map is strange and wont pickle onless it is at the top level
	lamda = params[0] # these params are per Wavelength, we could get more ganular parallelism with this!
	nk    = params[1]
	spectrum_list_generator=params[2].spectrum_list_generator
	parameter_list_generator=params[2].parameter_list_generator
	
	spectrum_list = spectrum_list_generator(lamda)
	sqrt_point_multiplicity = sqrt(len(spectrum_list))
	list_of_parameters = parameter_list_generator(lamda)
	#release the basic_setup file
	sub_error_list=[]
	if len(spectrum_list) > 0:
		for spectrum, parameters in zip (spectrum_list, list_of_parameters ):
			spectrum_calculated = TMM_spectrum_wrapper(nk_fit = nk, get_data=True, **parameters) # these parameters are per measurement, set of parameters for a single point tmm model
			error = (spectrum_calculated - spectrum)/sqrt_point_multiplicity # in the formulation this gets squared later, and we look at the per Wavelength rms error and normalizing it gives portability to weights
			#we'll put wieghting somewhere else so that the formulation is clear
			sub_error_list.append(error)
	else:
		print ("WARNING: No spectra at lamda = %f point, nk values here are only interpolations!" % lamda)
	

	return sub_error_list


def pointwise_rms_error_sum_wrapper(params): # for each Wavelength, wraps the previous function TRA_lamda_error to quikcly compute error spectrum
		from numpy import  sqrt
		point_error_list = spectrum_lamda_error(params)
		sum_err_square = 0.0
		for err in point_error_list:   #[0].[-1] will be error data of normal R and T
			sum_err_square += err**2
		lamda=params[0]
		spectrum_list_generator=params[2].spectrum_list_generator
		spectrum_list = spectrum_list_generator(lamda)
		sqrt_point_multiplicity = sqrt(len(spectrum_list))
		R_normal_error=point_error_list[0]*sqrt_point_multiplicity
		T_normal_error=point_error_list[-1]*sqrt_point_multiplicity
		return sqrt(sum_err_square),R_normal_error,T_normal_error



def rms_error_spectrum(lamda_list, nk_f, spectrum_list_generator, parameter_list_generator, threads = 0,
					input_data=None, test_setup=None, get_data=False):
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), test_setup ) )

	from multiprocessing import Pool, cpu_count
	if threads <= 0:
		threads = cpu_count()
	my_pool = Pool(threads)
	#print ('Using %i Threads' % threads)

	error_spectrum = my_pool.map(pointwise_rms_error_sum_wrapper, muh_inputs)
	#error_spectrum = map(pointwise_rms_error_sum_wrapper, muh_inputs)
	my_pool.terminate()
	if get_data:
		return error_spectrum
	else:
		temp_error=[]
		for tu in error_spectrum:
			temp_error.append(tu[0])
		error_spectrum=temp_error
		return error_spectrum

def sqr_rms_gradient_at_lamda(params, h_nk = 1e-6):
	from numpy import array
	def shift_nk(params, shift): # i use this code to make it clearer
		return [params[0], params[1]+shift, params[2]]

	e_p_0 = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk ))[0]**2
	e_m_0 = pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk ))[0]**2

	e_0_p = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*1.0j ))[0]**2
	e_0_m = pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk*1.0j ))[0]**2

	e_n = (e_p_0 - e_m_0)/(2.0*h_nk)
	e_k = (e_0_p - e_0_m)/(2.0*h_nk)
	gradient = array([e_n, e_k])
	return gradient


def sqr_rms_hessian_at_lamda(params, h_nk = 1e-4):
	from numpy import array
	def shift_nk(params, shift): # i use this code to make it clearer
		return [params[0], params[1]+shift, params[2]]
	e_0_0 = pointwise_rms_error_sum_wrapper(shift_nk(params,  0.0 ))[0]**2

	e_p_0 = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0+0.0j) ))[0]**2
	e_m_0 = pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk*( 1.0+0.0j) ))[0]**2

	e_0_p = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 0.0+1.0j) ))[0]**2
	e_0_m = pointwise_rms_error_sum_wrapper(shift_nk(params, -h_nk*( 0.0+1.0j) ))[0]**2

	e_p_p = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0+1.0j) ))[0]**2
	e_m_m = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0-1.0j) ))[0]**2

	e_p_m = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0-1.0j) ))[0]**2
	e_m_p = pointwise_rms_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0+1.0j) ))[0]**2

	e_nn = (e_p_0 - 2.0*e_0_0 + e_m_0)/(h_nk**2.0)
	e_kk = (e_0_p - 2.0*e_0_0 + e_0_m)/(h_nk**2.0)

	e_nk = (e_p_p - e_p_0 - e_0_p + 2.0*e_0_0 - e_m_0 - e_0_m + e_m_m)/(2.0*h_nk*h_nk) # two difference methods..., the first is supposed to be cheaper
	e_nk = (e_p_p - e_p_m - e_m_p + e_m_m)/(4.0*h_nk*h_nk)

	hessian = array([[e_nn, e_nk],[e_nk, e_kk]])
	return hessian

def pointwise_reducible_rms_error_sum_wrapper(params): # for each Wavelength, wraps the previous function TRA_lamda_error to quikcly compute error spectrum
	## calculates numeric derivatives
	#does some hessians and matrix match
	#from wikipedia

	
	e_0_0=pointwise_rms_error_sum_wrapper(params = params)[0]**2

	grad = sqr_rms_gradient_at_lamda(params = params)# h_nk = h_nk)
	hessian = sqr_rms_hessian_at_lamda(params = params)#, h_nk = h_nk)

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

def reducible_rms_error_spectrum(lamda_list, nk_f, spectrum_list_generator, parameter_list_generator, threads = 0,
								input_data=None, test_setup=None):
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), test_setup ) )

	from numpy import array, sqrt
	from multiprocessing import Pool, cpu_count
	if threads <= 0:
		threads = cpu_count()
	my_pool = Pool(threads)
	#print ('Using %i Threads' % threads)

	error_spectra = array(my_pool.map(pointwise_reducible_rms_error_sum_wrapper, muh_inputs)).T
	#error_spectra = array(map(pointwise_reducible_rms_error_sum_wrapper, muh_inputs)).T
	my_pool.terminate()

	reducible_error_spectrum = sqrt(error_spectra[0])
	irreducible_error_spectrum = sqrt(error_spectra[1])
	return reducible_error_spectrum, irreducible_error_spectrum



def single_lamda_rms_error_map(lamda, nlist, klist, spectrum_list_generator, parameter_list_generator, threads = 0):
	muh_inputs = []
	for n in nlist:
		for k in klist:
			muh_inputs.append( (lamda, n+1.0j*k, spectrum_list_generator, parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	if threads <= 0:
		threads = cpu_count()
	my_pool = Pool(threads)
	#print ('Using %i Threads' % threads)

	from numpy import reshape, array
	error_list = my_pool.map(pointwise_rms_error_sum_wrapper, muh_inputs)
	error_map = reshape( array(error_list), (len(nlist), len(klist) ))
	my_pool.terminate()
	return error_map #[nindex,kindex]


def spectrum_TMM_lamda(params): # This has to be at the top level because map is strange and wont pickle onless it is at the top level
	'''Returns a list of the spectrum values for each single point at a Wavelength'''
	lamda = params[0]
	nk    = params[1]
	parameter_list_generator = params[2]

	list_of_parameters = parameter_list_generator(lamda)
	spectrum_list = []
	for parameters in list_of_parameters :
		spectrum_list.append(  TMM_spectrum_wrapper(nk_fit = nk,  **parameters) )
	return spectrum_list


def TMM_spectra(lamda_list, nk_f,  parameter_list_generator, threads = 0):
	''''''
	#returns data in block like spectra[lamda][ spectrum] there are computational reasons why this order is this way
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	if threads <= 0:
		threads = cpu_count()
	my_pool = Pool(threads)
	#print ('Using %i Threads' % threads))

	spectra = my_pool.map(spectrum_TMM_lamda, muh_inputs)
	my_pool.terminate()
	#
	return spectra




###################
def fit_spectra_nk_sqr(lamda_list, spectrum_list_generator, parameter_list_generator,  
					nk_f_guess, delta_weight = 0.1, k_weight_fraction = 1.0, tolerance = 1e-5, 
					no_negative = True, interpolation_type = 'cubic', method = 'least_squares', threads = 0,
					input_data=None, test_setup=None, TB1=None
					):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi,exp,abs,sqrt, array, matmul, loadtxt, zeros, savetxt, inf, diff, ones
	from scipy.optimize import root, least_squares, minimize
	from TRANK import extrap


	#point_multiplicity = len(spectrum_list_generator(lamda_list[0]))
	#print(point_multiplicity)

	#point_multiplicity_list = [len(spectrum_list_generator(lamda)) for lamda in lamda_list ]
	#point_multiplicity = point_multiplicity_list[0]

	abs_delta_weight = sqrt(delta_weight**2  * (len(lamda_list)/(len(lamda_list)-1.0)))
    
	from multiprocessing import Pool, cpu_count
	if threads <= 0:
		threads = cpu_count()
	

	TB1.append ('Using %i Processes' % threads)
	
	my_pool=Pool(threads)
	
	
    
	def F_error(nk_list):
		
		c_nk_list = []
		muh_inputs = []
		for i in range(len(lamda_list)):
			nk = nk_list[2*i] + 1.0j*nk_list[2*i+1]				
			c_nk_list.append(nk)
			muh_inputs.append( (lamda_list[i], nk, test_setup) )
			'''
			ctx=mp.get_context('spawn')
			p=ctx.Process(target=spectrum_lamda_error,args=(child_conn,muh_inputs[i],))
			p.start()
			print (parent_conn.recv())
			p.join()
			
			'''
		
		#lock=RLock()	
		
		error_list_lists =my_pool.map(spectrum_lamda_error, muh_inputs)
		
		
		#error_list_lists = map(spectrum_lamda_error, muh_inputs)

		#combine the sub error lists into       #+=��mean to do so?
		error_list = []
		for sub_error_list in error_list_lists:
			error_list = error_list + sub_error_list
		
		
		delta_array = diff(c_nk_list)*abs_delta_weight

		error_list = error_list + list(delta_array.real) + list(k_weight_fraction * delta_array.imag)    #make up some points?
		
		return error_list

	####### now for a guess list ##############

	nk_guess_list = []
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		if no_negative:
			n = nk.real
			k = nk.imag
			if n < 0.0 : n = 0.0
			if k < 0.0 : k = 0.0
			nk_guess_list.append(n)
			nk_guess_list.append(k)
		else:
			nk_guess_list.append(nk.real)
			nk_guess_list.append(nk.imag)


	######### test
	if False: print(F_error(nk_guess_list))
	############ nk guessing over, time for creating and minimizing error function
	if method == 'least_squares':
		inputs = dict(fun = F_error,
					x0 = nk_guess_list,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2)
		if no_negative:
			inputs.update(dict( bounds = [zeros(2*len(lamda_list)),inf*ones(2*len(lamda_list))] ))
		solution = least_squares(**inputs ).x

	elif method == 'L-BFGS-B':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = nk_guess_list,
				method = 'L-BFGS-B',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(0,inf)]))
		solution = minimize(**inputs ).x

	elif method == 'SLSQP':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = nk_guess_list,
				method = 'SLSQP',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(0,inf)]))
		solution = minimize(**inputs ).x

	elif method == 'TNC':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = nk_guess_list,
				method = 'TNC',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(0,inf)]))
		solution = minimize(**inputs ).x

	else:
		raise ValueError("Invalid minimization method!")

	
	my_pool.terminate()
	my_pool.close()
	
		

	nk_list=[]
	for i in range(len(lamda_list)):
		nk_list.append(solution[2*i] + 1.0j*solution[2*i+1]  )

	fit_nk_f = extrap(lamda_list, nk_list, kind = interpolation_type)


	return fit_nk_f




def fit_spectra_nk_sqr_KK_compliant(lamda_list, lamda_fine, spectrum_list_generator, parameter_list_generator,  nk_f_guess,
								delta_weight = 0.1, k_weight_fraction = 1.0, tolerance = 1e-5, no_negative = True, 
								interpolation_type = 'cubic', method = 'least_squares', threads = 0,
								input_data=None, test_setup=None, TB1=None

								):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi,exp,abs,sqrt, array, matmul, loadtxt, zeros, savetxt, inf, diff, ones, mean
	from scipy.optimize import root, least_squares, minimize
	from TRANK import extrap
	#from TRANK import parallel_DKKT_n_from_lamda_k as KKT
	from TRANK import parallel_DKKT_n_from_lamda_k_with_edge_corrections as KKT


	#point_multiplicity = len(TR_pair_list_generator(lamda_list[0]))
	#print(point_multiplicity)

	abs_delta_weight = sqrt(delta_weight**2  * (len(lamda_list)/(len(lamda_list)-1.0)))

	from multiprocessing import  Pool
	if threads <= 0:
		threads = cpu_count()
	my_pool = Pool(threads)
	print ('Using %i Processes' % threads)
	def F_error(k_and_p_list):
		# the last value is the principle value
		#FYI -> k = array(k_and_p_list[0:-1])
		k = k_and_p_list[0:-1]
		p = k_and_p_list[-1]
		n = p + KKT(lamda_list = lamda_list, k = k, compute_pool = my_pool )



		muh_inputs = []
		for i in range(len(lamda_list)):
			nk = n[i]+1.0j*k[i]# double check this works properly later

			muh_inputs.append( (lamda_list[i], nk, test_setup ) )
		#print (zip(lamda_list, c_nk_list))
		error_list_lists = my_pool.map(spectrum_lamda_error, muh_inputs)
		#error_list_lists =my_pool.map(lamda_error, zip(lamda_list, c_nk_list))
		#print (error_list_lists)

		error_list = []
		for sub_error_list in error_list_lists:
			error_list = error_list + sub_error_list


		error_list = error_list + list( abs_delta_weight*diff(n) )   + list( k_weight_fraction * abs_delta_weight * diff(k))

		return error_list

	####### now for a guess list ##############

	guess_k_and_p_list= []
	p = 0.0
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		if no_negative and nk.imag < 0.0:
			k = 0
		else:
			k = nk.imag

		guess_k_and_p_list.append( k)
		p+= nk.real
	# now we put p at the end
	p = p/len(lamda_list)    # this is a guess for the principle value
	print ('principle value guess:',p)
	guess_k_and_p_list.append(p)


	######### test
	if False: print(F_error(guess_k_and_p_list)) #use this to see if the TR_error works
	############ nk guessing over, time for creating and minimizing error function

	if method == 'least_squares':
		inputs = dict(fun = F_error,
					x0 =  guess_k_and_p_list,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2)
		if no_negative:
			inputs.update(dict( bounds = [len(lamda_list)*[0.0]+[-inf],len(lamda_list)*[inf]+[inf]] ))
		solution = least_squares(**inputs ).x

	elif method == 'SLSQP':
		inputs = dict(fun = lambda x: 0.5*sum(array(F_error(x))**2),
				x0 = guess_k_and_p_list,
				method = 'SLSQP',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = len(lamda_list)*[(0,inf)]  + [(-inf,inf)] ))
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



	k_and_p_list = solution
	k = k_and_p_list[0:-1]
	p = k_and_p_list[-1]
	n = p +  KKT(lamda_list = lamda_list, k = k, compute_pool = my_pool )
	print ('Final principle value:',p)


	fit_nk_f = extrap(lamda_list, n + 1.0j*k, kind = interpolation_type)


	my_pool.terminate()
	my_pool.close()

	return fit_nk_f



























def TR(nk_fit, lamda, snell_angle_front, layer_index_of_fit,  nk_f_list,  thickness_list, coherency_list, tm_polarization_fraction): #this is just a fancy wrapper for inc_tmm

	nk_list = [ nk_f(lamda) for nk_f in nk_f_list]
	nk_list[layer_index_of_fit] = nk_fit # overwrite input nk_f with the fit one. basically, it would make the code much uglier if I made an excption

	te_result = inc_tmm('s', nk_list, thickness_list, coherency_list, snell_angle_front, lamda)
	tm_result = inc_tmm('p', nk_list, thickness_list, coherency_list, snell_angle_front, lamda)

	T = tm_polarization_fraction * tm_result['T'] + (1.0-tm_polarization_fraction) * te_result['T']
	R = tm_polarization_fraction * tm_result['R'] + (1.0-tm_polarization_fraction) * te_result['R']

	return T, R



def TRA_lamda(params): # This has to be at the top level because map is strange and wont pickle onless it is at the top level
	'''Returns a list of TRA for each single point at a Wavelength'''
	lamda = params[0]
	nk    = params[1]
	parameter_list_generator = params[2]

	list_of_parameters = parameter_list_generator(lamda)
	TRA_list = []
	for parameters in list_of_parameters :
		T, R = TR(nk_fit = nk,  **parameters)
		A = 1.0 - T - R
		TRA_list.append([T, R, A])
	return TRA_list


def TRA_lamda_error(params): # This has to be at the top level because map is strange and wont pickle onless it is at the top level
	lamda = params[0] # these params are per Wavelength
	nk    = params[1]
	TR_pair_list_generator = params[2]
	parameter_list_generator = params[3]

	TR_pair_list = TR_pair_list_generator(lamda)
	list_of_parameters = parameter_list_generator(lamda)
	sub_error_list=[]
	for TR_pair, parameters in zip (TR_pair_list, list_of_parameters ):
		Tc,Rc = TR(nk_fit = nk,  **parameters) # these parameters are per measurement, set of parameters for a single point tmm model
		Ac = 1.0 - Tc - Rc
		A_pair = 1.0 - TR_pair[0] - TR_pair[1]
		sub_error_list.append(TR_pair[0]-Tc)
		sub_error_list.append(TR_pair[1]-Rc)
		sub_error_list.append(A_pair    -Ac)

	return sub_error_list







###################
def fit_TRA_nk_sqr(lamda_list, TR_pair_list_generator, parameter_list_generator,  nk_f_guess, delta_weight = 0.1, tolerance = 1e-4, no_negative = True, interpolation_type = 'cubic', method = 'least_squares'):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi,exp,abs,sqrt, array, matmul, loadtxt, zeros, savetxt, inf, diff, ones
	from scipy.optimize import root, least_squares, minimize
	from TRANK import extrap


	point_multiplicity = len(TR_pair_list_generator(lamda_list[0]))
	#print(point_multiplicity)
	# 3.0  is from T, R, and A
	abs_delta_weight = sqrt(delta_weight**2  * point_multiplicity * 3.0 * (len(lamda_list)/(len(lamda_list)-1.0)))

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	#my_pool = Pool()
	#my_pool = Pool(1)


	def TR_error(nk_list):

		c_nk_list = []
		muh_inputs = []
		for i in range(len(lamda_list)):
			nk = nk_list[2*i] + 1.0j*nk_list[2*i+1]
			c_nk_list.append(nk)
			muh_inputs.append( (lamda_list[i], nk, TR_pair_list_generator, parameter_list_generator ) )
			
		#print (zip(lamda_list, c_nk_list))
		error_list_lists = my_pool.map(TRA_lamda_error, muh_inputs)
		#error_list_lists =my_pool.map(lamda_error, zip(lamda_list, c_nk_list))
		#print (error_list_lists)

		error_list = []
		for sub_error_list in error_list_lists:
			error_list = error_list + sub_error_list

		delta_array = diff(c_nk_list)*abs_delta_weight
		#delta_errors =(delta_array.real**2 + delta_array.imag**2) * abs_delta_weight

		#error_list = error_list + list(delta_errors)
		error_list = error_list + list(delta_array.real) + list(delta_array.imag)

		return array(error_list)

	####### now for a guess list ##############

	nk_guess_list = []
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		nk_guess_list.append(abs(nk.real))
		nk_guess_list.append(abs(nk.imag))

	######### test
	#print(TR_error(nk_guess_list))
	############ nk guessing over, time for creating and minimizing error function
	if method == 'least_squares':
		inputs = dict(fun = TR_error,
					x0 = nk_guess_list,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2)
		if no_negative:
			inputs.update(dict( bounds = [zeros(2*len(lamda_list)),inf*ones(2*len(lamda_list))] ))
		solution = least_squares(**inputs ).x

	elif method == 'L-BFGS-B':
		inputs = dict(fun = lambda x: 0.5*sum(TR_error(x)**2),
				x0 = nk_guess_list,
				method = 'L-BFGS-B',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(0,inf)]))
		solution = minimize(**inputs ).x

	elif method == 'SLSQP':
		inputs = dict(fun = lambda x: 0.5*sum(TR_error(x)**2),
				x0 = nk_guess_list,
				method = 'SLSQP',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(0,inf)]))
		solution = minimize(**inputs ).x

	elif method == 'TNC':
		inputs = dict(fun = lambda x: 0.5*sum(TR_error(x)**2),
				x0 = nk_guess_list,
				method = 'TNC',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = 2*len(lamda_list)*[(0,inf)]))
		solution = minimize(**inputs ).x

	else:
		raise ValueError("Invalid minimization method!")

	my_pool.terminate()
	my_pool.close()
	n_list=[]
	k_list=[]
	for i in range(len(lamda_list)):
		n_list.append(solution[2*i]  )
		k_list.append(solution[2*i+1])

	nf = extrap(lamda_list, n_list, kind = interpolation_type)
	kf = extrap(lamda_list, k_list, kind = interpolation_type)

	def fit_nk_f(lamda):
		return nf(lamda) + 1.0j*kf(lamda)

	return fit_nk_f













def fit_TRA_nk_sqr_KK_compliant(lamda_list, lamda_fine, TR_pair_list_generator, parameter_list_generator,  nk_f_guess,
								delta_weight = 0.1, tolerance = 1e-5, no_negative = True, interpolation_type = 'cubic', method = 'least_squares'):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi,exp,abs,sqrt, array, matmul, loadtxt, zeros, savetxt, inf, diff, ones, mean
	from scipy.optimize import root, least_squares, minimize
	from TRANK import extrap


	point_multiplicity = len(TR_pair_list_generator(lamda_list[0]))
	#print(point_multiplicity)
	# 3.0  is from T, R, and A
	abs_delta_weight = sqrt(delta_weight**2  * point_multiplicity * 3.0 * (len(lamda_list)/(len(lamda_list)-1.0)))

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	#my_pool = Pool()
	#my_pool = Pool(1)


	def TR_error(k_and_p_list):
		# the last value is the principle value
		#FYI -> k = array(k_and_p_list[0:-1])
		k = k_and_p_list[0:-1]
		p = k_and_p_list[-1]
		n = p + KK_lamda(lamda_list = lamda_list, lamda_fine = lamda_fine,  k = k )



		muh_inputs = []
		for i in range(len(lamda_list)):
			nk = n[i]+1.0j*k[i]# double check this works properly later

			muh_inputs.append( (lamda_list[i], nk, TR_pair_list_generator, parameter_list_generator ) )

		#print (zip(lamda_list, c_nk_list))
		error_list_lists = my_pool.map(TRA_lamda_error, muh_inputs)
		#error_list_lists =my_pool.map(lamda_error, zip(lamda_list, c_nk_list))
		#print (error_list_lists)

		error_list = []
		for sub_error_list in error_list_lists:
			error_list = error_list + sub_error_list


		error_list = error_list + list( abs_delta_weight*diff(n) )   + list( abs_delta_weight * diff(k))

		return error_list

	####### now for a guess list ##############

	guess_k_and_p_list= []
	p = 0.0
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		if no_negative and nk.imag < 0.0:
			k = 0
		else:
			k = nk.imag

		guess_k_and_p_list.append( k)
		p+= nk.real
	# now we put p at the end
	p = p/len(lamda_list) - 1.0   # this is a guess for the principle value
	print ('principle value guess:',p)
	guess_k_and_p_list.append(p)


	######### test
	if False: print(TR_error(guess_k_and_p_list)) #use this to see if the TR_error works
	############ nk guessing over, time for creating and minimizing error function

	if method == 'least_squares':
		inputs = dict(fun = TR_error,
					x0 =  guess_k_and_p_list,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2)
		if no_negative:
			inputs.update(dict( bounds = [len(lamda_list)*[0.0]+[-inf],len(lamda_list)*[inf]+[inf]] ))
		solution = least_squares(**inputs ).x

	elif method == 'SLSQP':
		inputs = dict(fun = lambda x: 0.5*sum(TR_error(x)**2),
				x0 = guess_k_and_p_list,
				method = 'SLSQP',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = len(lamda_list)*[(0,inf)]  + [(-inf,inf)] ))
		solution = minimize(**inputs ).x

	elif method == 'L-BFGS-B':
		inputs = dict(fun = lambda x: 0.5*sum(TR_error(x)**2),
				x0 = guess_k_and_p_list,
				method = 'L-BFGS-B',
				tol = tolerance,
				options = {'disp' : True, 'iprint': 2}  )
		if no_negative:
			inputs.update(dict( bounds = len(lamda_list)*[(0,inf)]  + [(-inf,inf)] ))
		solution = minimize(**inputs ).x

	elif method == 'TNC':
		inputs = dict(fun = lambda x: 0.5*sum(TR_error(x)**2),
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

	nf = extrap(lamda_list, n, kind = interpolation_type)
	kf = extrap(lamda_list, k, kind = interpolation_type)

	def fit_nk_f(lamda):
		return nf(lamda) + 1.0j*kf(lamda)

	return fit_nk_f



#########################################
def fit_TRA_epsilon_sqr(lamda_list, TR_pair_list_generator, parameter_list_generator,  nk_f_guess, delta_weight = 0.1, tolerance = 1e-4, no_negative = True):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi,exp,abs,sqrt, array, matmul, loadtxt, zeros, savetxt, inf, diff, ones
	from scipy.optimize import root, least_squares


	point_multiplicity = len(TR_pair_list_generator(lamda_list[0]))
	#print(point_multiplicity)
	# 3.0  is from T, R, and A
	abs_delta_weight = sqrt(delta_weight**2  * point_multiplicity * 3.0 * (len(lamda_list)/(len(lamda_list)-1.0)))

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	#my_pool = Pool()
	#my_pool = Pool(1)


	def TR_error(nk_list):

		c_nk_list = zeros(len(lamda_list),dtype = complex)
		muh_inputs = []
		for i in range(len(lamda_list)):
			nk = nk_list[2*i] + 1.0j*nk_list[2*i+1]
			c_nk_list[i]=nk
			muh_inputs.append( (lamda_list[i], nk, TR_pair_list_generator, parameter_list_generator ) )

		#print (zip(lamda_list, c_nk_list))
		error_list_lists = my_pool.map(TRA_lamda_error, muh_inputs)
		#error_list_lists =my_pool.map(lamda_error, zip(lamda_list, c_nk_list))
		#print (error_list_lists)

		error_list = []
		for sub_error_list in error_list_lists:
			error_list = error_list + sub_error_list

		#delta_array = diff(c_nk_list)*abs_delta_weight

		delta_array = diff(c_nk_list**2)*abs_delta_weight
		#delta_errors =(delta_array.real**2 + delta_array.imag**2) * abs_delta_weight

		#error_list = error_list + list(delta_errors)
		error_list = error_list + list(delta_array.real) + list(delta_array.imag)

		return error_list

	####### now for a guess list ##############

	nk_guess_list = []
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		nk_guess_list.append(abs(nk.real))
		nk_guess_list.append(abs(nk.imag))

	######### test
	#print(TR_error(nk_guess_list))
	############ nk guessing over, time for creating and minimizing error function
	if no_negative:
		solution = least_squares(TR_error,
					x0 = nk_guess_list,
					bounds = [zeros(2*len(lamda_list)),inf*ones(2*len(lamda_list))] ,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2).x
	else:
		solution = least_squares(TR_error,
					x0 = nk_guess_list,
					ftol = tolerance,
					xtol = tolerance,
					gtol = tolerance,
					verbose = 2).x

	my_pool.terminate()
	my_pool.close()
	n_list=[]
	k_list=[]
	for i in range(len(lamda_list)):
		n_list.append(solution[2*i]  )
		k_list.append(solution[2*i+1])

	nf = extrap(lamda_list, n_list, kind = 'cubic')
	kf = extrap(lamda_list, k_list, kind = 'cubic')

	def fit_nk_f(lamda):
		return nf(lamda) + 1.0j*kf(lamda)

	return fit_nk_f




###################
def fit_TRA_holy_nksqr_BFGS(lamda_list, TR_pair_list_generator, parameter_list_generator,  nk_f_guess, delta_weight = 0.1, tolerance = 1e-4, no_negative = True):
	'''n_front and n_back must be real valued for this to work without caveats.
thickness and lambda can be any units, so long as they are the same, lamda_list must be sorted'''


	from numpy import pi,exp,abs,sqrt, array, matmul, loadtxt, zeros, savetxt, inf, diff, ones
	from scipy.optimize import root, least_squares, minimize
	from TRANK import extrap


	point_multiplicity = len(TR_pair_list_generator(lamda_list[0]))
	#print(point_multiplicity)
	# 3.0  is from T, R, and A
	abs_delta_weight = sqrt(delta_weight**2  * point_multiplicity * 3.0 * (len(lamda_list)/(len(lamda_list)-1.0)))

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())


	def TR_error(nk_list):

		c_nk_list = []
		muh_inputs = []
		for i in range(len(lamda_list)):
			nk = nk_list[2*i] + 1.0j*nk_list[2*i+1]
			c_nk_list.append(nk)
			muh_inputs.append( (lamda_list[i], nk, TR_pair_list_generator, parameter_list_generator ) )

		#print (zip(lamda_list, c_nk_list))
		error_list_lists = my_pool.map(TRA_lamda_error, muh_inputs)
		#error_list_lists =my_pool.map(lamda_error, zip(lamda_list, c_nk_list))
		#print (error_list_lists)

		#error_list = []
		#for sub_error_list in error_list_lists:
		#	error_list = error_list + sub_error_list

		base_line_error_square = (array(error_list_lists)**2).sum()


		delta_array = diff(c_nk_list)*abs_delta_weight
		#delta_errors =(delta_array.real**2 + delta_array.imag**2) * abs_delta_weight

		#error_list = error_list + list(delta_errors)
		#error_list = error_list + list(delta_array.real) + list(delta_array.imag)

		delta_errors_square = (delta_array.real**2 + delta_array.imag**2).sum() * abs_delta_weight

		return base_line_error_square + delta_errors_square

	####### now for a guess list ##############

	nk_guess_list = []
	for i in range(len(lamda_list)):
		nk = nk_f_guess(lamda_list[i])
		nk_guess_list.append(abs(nk.real))
		nk_guess_list.append(abs(nk.imag))

	######### test
	#print(TR_error(nk_guess_list))
	############ nk guessing over, time for creating and minimizing error function
	if no_negative:
		solution = minimize(TR_error,
					x0 = nk_guess_list,
					method = 'L-BFGS-B',
					bounds = 2*len(lamda_list)*[(0,inf)] ,
					tol = tolerance,
					options = {'disp' : True } ).x
	else:
		solution = minimize(TR_error,
					x0 = nk_guess_list,
					method = 'L-BFGS-B',
					tol = tolerance,
					options = {'disp' : True } ).x

	my_pool.terminate()
	my_pool.close()
	n_list=[]
	k_list=[]
	for i in range(len(lamda_list)):
		n_list.append(solution[2*i]  )
		k_list.append(solution[2*i+1])

	nf = extrap(lamda_list, n_list, kind = 'cubic')
	kf = extrap(lamda_list, k_list, kind = 'cubic')

	def fit_nk_f(lamda):
		return nf(lamda) + 1.0j*kf(lamda)

	return fit_nk_f









def TRA_spectra(lamda_list, nk_f,  parameter_list_generator):
	''''''
	#returns data in block like spectra[lamda][ param][ spectrum] there are computational reasons why this order is this way
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	spectra = my_pool.map(TRA_lamda, muh_inputs)
	my_pool.terminate()
	#
	return spectra




def pointwise_TRA_error_sum_wrapper(params): # for each Wavelength, wraps the previous function TRA_lamda_error to quikcly compute error spectrum
		from numpy import  sqrt
		point_error_list = TRA_lamda_error(params)
		sum_err_square = 0.0
		for err in point_error_list:
			sum_err_square += err**2
		return sqrt(sum_err_square/len(point_error_list))



def rms_TRA_error_spectrum(lamda_list, nk_f, TR_pair_list_generator, parameter_list_generator):
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), TR_pair_list_generator, parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	error_spectrum = my_pool.map(pointwise_TRA_error_sum_wrapper, muh_inputs)

	my_pool.terminate()
	return error_spectrum

def gradient_at_lamda(params, h_nk = 1e-6):
	from numpy import array
	def shift_nk(params, shift): # i use this code to make it clearer
		return [params[0], params[1]+shift, params[2], params[3]]

	e_p_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk ))**2
	e_m_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params, -h_nk ))**2

	e_0_p = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*1.0j ))**2
	e_0_m = pointwise_TRA_error_sum_wrapper(shift_nk(params, -h_nk*1.0j ))**2

	e_n = (e_p_0 - e_m_0)/(2.0*h_nk)
	e_k = (e_0_p - e_0_m)/(2.0*h_nk)
	gradient = array([e_n, e_k])
	return gradient


def hessian_at_lamda(params, h_nk = 1e-4):
	from numpy import array
	def shift_nk(params, shift): # i use this code to make it clearer
		return [params[0], params[1]+shift, params[2], params[3]]
	e_0_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params,  0.0 ))**2

	e_p_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0+0.0j) ))**2
	e_m_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params, -h_nk*( 1.0+0.0j) ))**2

	e_0_p = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*( 0.0+1.0j) ))**2
	e_0_m = pointwise_TRA_error_sum_wrapper(shift_nk(params, -h_nk*( 0.0+1.0j) ))**2

	e_p_p = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0+1.0j) ))**2
	e_m_m = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0-1.0j) ))**2

	e_p_m = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0-1.0j) ))**2
	e_m_p = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0+1.0j) ))**2

	e_nn = (e_p_0 - 2.0*e_0_0 + e_m_0)/(h_nk**2.0)
	e_kk = (e_0_p - 2.0*e_0_0 + e_0_m)/(h_nk**2.0)

	e_nk = (e_p_p - e_p_0 - e_0_p + 2.0*e_0_0 - e_m_0 - e_0_m + e_m_m)/(2.0*h_nk*h_nk) # two difference methods..., the first is supposed to be cheaper
	e_nk = (e_p_p - e_p_m - e_m_p + e_m_m)/(4.0*h_nk*h_nk)

	hessian = array([[e_nn, e_nk],[e_nk, e_kk]])
	return hessian

def pointwise_reducible_TRA_error_sum_wrapper(params): # for each Wavelength, wraps the previous function TRA_lamda_error to quikcly compute error spectrum
	## calculates numeric derivatives
	#does some hessians and matrix match
	#from wikipedia

	if False: #doing all the stencils in one go is going to be much faster than call the individual functions, but this is not as easy to debug
		def shift_nk(params, shift): # i use this code to make it clearer
			return [params[0], params[1]+shift, params[2], params[3]]

		# building stencil components
		e_0_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params,  0.0 ))**2
		e_p_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk ))**2
		e_m_0 = pointwise_TRA_error_sum_wrapper(shift_nk(params, -h_nk ))**2

		e_0_p = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*1.0j ))**2
		e_0_m = pointwise_TRA_error_sum_wrapper(shift_nk(params, -h_nk*1.0j ))**2

		e_p_p = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0+1.0j) ))**2
		e_m_m = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0-1.0j) ))**2

		e_p_m = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*( 1.0-1.0j) ))**2
		e_m_p = pointwise_TRA_error_sum_wrapper(shift_nk(params,  h_nk*(-1.0+1.0j) ))**2

		# building derivatives
		e_n = (e_p_0 - e_m_0)/(2.0*h_nk)
		e_k = (e_0_p - e_0_m)/(2.0*h_nk)

		e_nn = (e_p_0 - 2*e_0_0 + e_m_0)/(h_nk**2)
		e_kk = (e_0_p - 2*e_0_0 + e_0_m)/(h_nk**2)

		e_nk = (e_p_p - e_p_0 - e_0_p + 2.0*e_0_0 - e_m_0 - e_0_m + e_m_m)/(2.0*h_nk*h_nk) # two difference methods..., the first is supposed to be cheaper
		e_nk = (e_p_p - e_p_m - e_m_p + e_m_m)/(4.0*h_nk*h_nk)

		# now for the gradient and hessian
		grad = array([e_n, e_k])
		hessian = array([[e_nn, e_nk],[e_nk, e_kk]])
	else:
		e_0_0 = pointwise_TRA_error_sum_wrapper(params = params)**2
		grad = gradient_at_lamda(params = params)# h_nk = h_nk)
		hessian = hessian_at_lamda(params = params)#, h_nk = h_nk)

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

def reducible_rms_TRA_error_spectrum(lamda_list, nk_f, TR_pair_list_generator, parameter_list_generator):
	muh_inputs = []
	for lamda in lamda_list:
		muh_inputs.append( (lamda, nk_f(lamda), TR_pair_list_generator, parameter_list_generator ) )

	from numpy import array, sqrt
	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())
	error_spectra = array(my_pool.map(pointwise_reducible_TRA_error_sum_wrapper, muh_inputs)).T
	my_pool.terminate()

	reducible_error_spectrum = sqrt(error_spectra[0])
	irreducible_error_spectrum = sqrt(error_spectra[1])
	return reducible_error_spectrum, irreducible_error_spectrum



def single_lamda_TRA_error_map(lamda, nlist, klist, TR_pair_list_generator, parameter_list_generator):
	muh_inputs = []
	for n in nlist:
		for k in klist:
			muh_inputs.append( (lamda, n+1.0j*k, TR_pair_list_generator, parameter_list_generator ) )

	from multiprocessing import Pool, cpu_count
	my_pool = Pool(cpu_count())

	from numpy import reshape, array
	error_list = my_pool.map(pointwise_TRA_error_sum_wrapper, muh_inputs)
	error_map = reshape( array(error_list), (len(nlist), len(klist) ))
	my_pool.terminate()
	return error_map #[nindex,kindex]

def find_min_indices_2d_array(thing):
	min_indices = [0,0]
	for i in range(thing.shape[0]):
		for j in range(thing.shape[1]):
			if thing[i][j] < thing[min_indices[0]][min_indices[1]]:
				min_indices = [i,j]
	return min_indices










###### end lib
