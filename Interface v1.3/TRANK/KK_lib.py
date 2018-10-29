



def dual_grid_direct_KK_n_from_lamda_k(lamda_list, lamda_fine,  k,  cshift = 1e-4) :
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





######### DKKT
from numpy import abs, pi
from numpy import log as ln
def g(x,y):
	cshift = 0.000001j
	return (x+y)*ln(abs(x+y)) + (x-y)*ln(abs(x-y)+cshift)


def n_jj_omega(jj, omega, omega_list): # the index kernel
	term1 = g(omega,omega_list[jj-1])/(omega_list[jj] - omega_list[jj-1])
	term2 = - ( (omega_list[jj+1] - omega_list[jj-1]) * g(omega,omega_list[jj]) )/( (omega_list[jj] - omega_list[jj-1]) * (omega_list[jj+1] - omega_list[jj]))
	term3 = g(omega,omega_list[jj+1])/(omega_list[jj+1] - omega_list[jj])
	return (term1 + term2 + term3).real/pi



def DKKT_n_from_lamda_k(lamda_list, k):
	n_list = zeros(len(lamda_list))
	omega_list = 1.0/array(lamda_list)
	for i in range(len(lamda_list)):
		lamda = lamda_list[i]

		n = 0
		for jj in range(1,len(lamda_list)-1):
			n += n_jj_omega(jj = jj, omega = 1.0/lamda_list[i], omega_list = omega_list  ) * k[jj]

		n_list[i] = n

	return n_list





def single_point_DKKT_from_omega_k(inputs):
	point_index = inputs['point_index']
	omega_list = inputs['omega_list']
	k = inputs['k']

	n = 0
	for jj in range(1,len(omega_list)-1):
		n += n_jj_omega(jj = jj, omega = omega_list[point_index], omega_list = omega_list  ) * k[jj]
	return n

from numpy import array
def parallel_DKKT_n_from_lamda_k(lamda_list, k, compute_pool):
	# we could build the KK transform matrix out of the lambda points and
	# not have to recompute it later. but this is very clear for now

	omega_list = 1.0/array(lamda_list)
	input_list = []
	for i in range(len(lamda_list)):
		input_list.append({'point_index': i, 'omega_list': omega_list, 'k': k})

	n_list = compute_pool.map(single_point_DKKT_from_omega_k, input_list )
	return n_list



def upper_bound_extrapolation_order_0(lamda_list, k):
	'''Requires k to be lambda sorted'''
	lamda_list = array(lamda_list)
	k_b = k[-1]
	lamda_b = lamda_list[-1]

	#delta = (lamda_list[-1] - lamda_list[-2])
	delta = 10.0
	eta = 0.01
	from numpy import log as ln
	from numpy import pi
	#print (lamda_b, k_b)
	n_correction =  k_b/(pi*(1-eta)) *ln( abs(  (lamda_b+delta)**2  - (1-eta)* lamda_list**2 )/(lamda_b )**2 )

	return (n_correction.real)


def upper_bound_extrapolation_order_1_old(lamda_list, k):
	'''Requires k to be lambda sorted'''
	#k_b = k[-1]
	k_b = 0.0 #really bad about the edge, so we disable it and only use the derivative
	lamda_b = lamda_list[-1]
	k_prime =  (k[-1] - k[-2])/(lamda_list[-1] - lamda_list[-2])

	lamda_list = array(lamda_list)
	#print (lamda_b, k_b, k_prime)
	from numpy import log as ln
	from numpy import pi

	part_1 = (k_b + k_prime * (lamda_list - lamda_b) ) * ln( abs(lamda_list - lamda_b + 0.0000001j) )
	part_2 = (k_b - k_prime * (lamda_list + lamda_b) ) * ln( abs(lamda_list + lamda_b) )
	part_3 = 2 * ( k_prime*lamda_b - k_b ) * ln( abs( lamda_b ) )

	n_correction =  1.0/pi *(part_1 + part_2 + part_3)

	return (n_correction.real)


def upper_bound_extrapolation_order_1(lamda_list, k):
	'''Requires k to be lambda sorted'''
	#k_b = k[-1]
	lamda_b = lamda_list[-1]
	k_prime =  (k[-1] - k[-2])/(lamda_list[-1] - lamda_list[-2])

	lamda_list = array(lamda_list)
	#print (lamda_b, k_b, k_prime)
	from numpy import log as ln
	from numpy import pi

	part_1 =    (lamda_list - lamda_b) * ln( abs(lamda_list - lamda_b + 0.0000001j) ) # +0.0000001j
	part_2 =  - (lamda_list + lamda_b) * ln( abs(lamda_list + lamda_b) )
	part_3 = 2 * lamda_b * ln( abs( lamda_b ) )

	n_correction =  - k_prime / pi *(part_1 + part_2 + part_3)

	return (n_correction.real)




def parallel_DKKT_n_from_lamda_k_with_edge_corrections(lamda_list, k, compute_pool):
	# we could build the KK transform matrix out of the lambda points and
	# not have to recompute it later. but this is very clear for now

	n_list = array(parallel_DKKT_n_from_lamda_k(lamda_list, k, compute_pool))

	ub1 = upper_bound_extrapolation_order_1(lamda_list, k)
	ub0 = upper_bound_extrapolation_order_0(lamda_list, k)
	return n_list + 1.0*ub0 + 0.5*ub1
