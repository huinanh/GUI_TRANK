


def compute_outer_log_levels(error_map, base10delta = 0.2):
	from numpy import floor, ceil, log10, logspace
	def bottom_log(x):
		return base10delta * floor( log10(x)/base10delta )

	def top_log(x):
		return base10delta * ceil( log10(x)/base10delta )

	low =  bottom_log( error_map.min() )
	high = top_log(    error_map.max() )

	nlevels = (high-low)/base10delta +1
	levels =  logspace(low, high ,nlevels  )

	return levels


def compute_inner_log_levels(error_map, base10delta = 1.0):
	from numpy import floor, ceil, log10, logspace
	def bottom_log(x):
		return base10delta * ceil( log10(x)/base10delta )

	def top_log(x):
		return base10delta * floor( log10(x)/base10delta )

	low =  bottom_log( error_map.min() )
	high = top_log(    error_map.max() )

	nlevels = (high-low)/base10delta +1
	levels =  logspace(low, high ,nlevels  )

	return levels





from TRANK import single_lamda_rms_error_map, functionize_nk_file, try_mkdir, find_min_indices_2d_array, extrap, rms_error_spectrum
from numpy import ceil, floor, loadtxt

if __name__=='__main__':

	data_directory = 'TRANK_nk_fit/'
	map_direct = 'TRANK_error_maps/'

	from basic_setup  import fit_nk_f, spectrum_list_generator,   parameter_list_generator
	fit_nk_f =  functionize_nk_file( data_directory+'fit_nk_fine.txt', skiprows = 0)
	lamda_list = loadtxt(data_directory+'fit_nk.txt' , unpack = True, usecols = [0])
	lamda_fine = loadtxt(data_directory+'fit_nk_fine.txt' , unpack = True, usecols = [0])

	rms_error_fine = rms_error_spectrum(lamda_list = lamda_fine, nk_f = fit_nk_f,
										spectrum_list_generator = spectrum_list_generator,
										parameter_list_generator = parameter_list_generator)
	rms_f = extrap(lamda_fine, rms_error_fine, kind = 'linear')



	try_mkdir(map_direct )

	from matplotlib.backends.backend_pdf import PdfPages
	pdf_page_combined = PdfPages(map_direct + 'all_error_maps.pdf')


	from matplotlib.pylab import *
	from matplotlib.colors import LogNorm
	from matplotlib import ticker

	fit_n_max = fit_nk_f(lamda_list).real.max()
	fit_k_max = fit_nk_f(lamda_list).imag.max()

	num_n_points = 50
	num_k_points = 100

	nmin, nmax = 0.01,  fit_n_max *2.0
	kmin, kmax = 0.0,  fit_k_max *1.1

	#dn = 0.01
	#dk = 0.01
	use_logscale = True

	nlist = linspace(nmin, nmax, num_n_points)
	klist = linspace(kmin, kmax, num_k_points)


	coarse_lamda_list = lamda_list # uses the spacing of mesh points
	#coarse_lamda_list = arange(min(lamda_fine),max(lamda_fine)+.0001, 50) # fixed spacing
	#coarse_lamda_list = [600] # single point
	for lamda in coarse_lamda_list:
		print('lambda :',lamda)




		error_map = single_lamda_rms_error_map(lamda = lamda,
				nlist = nlist,
				klist = klist,
				spectrum_list_generator = spectrum_list_generator,
				parameter_list_generator = parameter_list_generator) *100.0





		f = figure(figsize =(3.2,2.5), dpi = 220*2.0/3.0 * 2.0)
		N, K = meshgrid(nlist, klist, indexing = 'ij')

		ylabel('k', fontsize = 10)
		xlabel('n',fontsize = 10)
		title('$\lambda$ = %i nm' % lamda)
		minorticks_on()
		gca().tick_params(axis='both', which='major', labelsize=10)

		if use_logscale:

			cf_levels  = compute_outer_log_levels(error_map, base10delta = 0.1)
			CF = contourf(N,K, error_map,  norm = LogNorm(), levels = cf_levels , locator=ticker.LogLocator() )
			cbar = colorbar(CF, ticks = compute_inner_log_levels(error_map, base10delta = 0.2))
			cbar.ax.set_ylabel('RMS Error (%)')
			#cbar.ax.yaxis.set_ticks(log10(cf_levels) , minor=True)
		else:
			cf_level_delta = 0.5
			error_min = error_map.min()
			error_max = error_map.max()
			cfmax = (error_max - error_min) *0.15 +error_min
			cf_levels = arange( floor(error_min), ceil(cfmax) ,cf_level_delta  )

			CF = contourf(N,K, error_map, levels = cf_levels )
			cbar = colorbar(CF)
			cbar.ax.set_ylabel('RMS Error (%)')

		min_indices = find_min_indices_2d_array(error_map)
		global_min_rms = error_map[min_indices[0],min_indices[1]]
		plot(nlist[min_indices[0]], klist[min_indices[1]], marker = 'o', markerfacecolor = 'white',markeredgecolor = 'white',linewidth = 0.0, markersize = 2,  label = 'On Map %.2e %%'%global_min_rms)

		rms_fit = rms_f(lamda) *100
		plot(fit_nk_f(lamda).real, fit_nk_f(lamda).imag, marker = '+', markerfacecolor = 'grey',markeredgecolor = 'grey',linewidth = 0.0,  label = 'Fit %.2e %%'%rms_fit)


		legend(fontsize = 6, handlelength = 1, handletextpad = 0.2)



		gcf().tight_layout(pad=0.15)
		pdf_page_combined.savefig(  transparent = True)
		close(f)

	pdf_page_combined.close()
