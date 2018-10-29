


from .TRANK_lib import (fit_spectra_nk_sqr,
						fit_spectra_nk_sqr_KK_compliant,
						rms_error_spectrum,
						reducible_rms_error_spectrum,
						sqr_rms_gradient_at_lamda,
						sqr_rms_hessian_at_lamda,
						pointwise_TRA_error_sum_wrapper,
						TMM_spectrum_wrapper,
						TMM_spectra,
						single_lamda_rms_error_map,
						find_min_indices_2d_array )

from .NKT_lib import ( NKT_fit_spectra_nk_sqr, NKT_fit_spectra_nk_sqr_KK_compliant, NKT_rms_error_spectrum, NKT_reducible_rms_error_spectrum)
from .NKT_iterator import ( NKT_error_adaptive_iterative_fit_spectra)



from .KK_lib import (dual_grid_direct_KK_n_from_lamda_k, DKKT_n_from_lamda_k,
					parallel_DKKT_n_from_lamda_k,
					upper_bound_extrapolation_order_0, upper_bound_extrapolation_order_1, parallel_DKKT_n_from_lamda_k_with_edge_corrections)


from .helpers import (extrap, extrap_c, functionize_nk_file, functionize_frequency_and_permittivity_file, try_mkdir, nk_plot, error_plot, compute_coarse_and_fine_grid)

from .iterator import error_adaptive_iterative_fit
from .iterator import error_adaptive_iterative_fit_spectra




from .TRANK_lib import (fit_TRA_nk_sqr,
						fit_TRA_nk_sqr_KK_compliant,
						rms_TRA_error_spectrum,
						reducible_rms_TRA_error_spectrum,
						gradient_at_lamda,
						hessian_at_lamda,
						pointwise_TRA_error_sum_wrapper,
						TR,
						TRA_spectra,
						single_lamda_TRA_error_map,
						find_min_indices_2d_array )
