#Auxilary functions
import curvefunc
"""
FUNCTIONS
    calc_alive_time(curve, limit)
    get_lightcurve_ULX_flux_limit(curve_0, Lx)
    load_curve_file(path)
    load_curve_file_skip_footer(path)
    multiply_curve(curve, period, time_range)
    scale_light_curve_period(curve, original_period, new_period)
"""

import batchfunc
"""
FUNCTIONS
    make_xcm(xcm_filename, parameters, lightcurve_filename)
    run_ulxlc(xcm_filename, parameters, lightcurve_filename)
    run_xcm(xcm_filename)
"""


import process_startack

"""
FUNCTIONS
    main()
"""



import process_startrack_new

import ulxlc_runner
"""
FUNCTIONS
    main()
    wrapper(id_dincl_tuple)
"""

import create_normalisation_lookup
"""
FUNCTIONS
    calc_N_lim(c0_max, Lx)
    main()
"""

import curve_classifier
"""
FUNCTIONS
    classify(N_lim, lc_min, lc_max)
    main()
"""

import create_erass_curves
"""
FUNCTIONS
    create_parent_population()
    get_curve(row)
    main()
"""

import erass_simulation
"""
FUNCTIONS
    erass_sample_curve(curve, curve_period, N_lim, number_of_repeats)
    get_curve(row)
    main()
"""


import process_erass_curves
"""
FUNCTIONS
    erass_sample_curve(curve, curve_period, N_lim, number_of_repeats)
    filename_from_row(row)
    main()
"""





# Legacy
import alive_dead_calc
import bh_vs_classification
import eROSITA_lightcurve_sampler
import eROSITA_results_processor

#other
import appendix_tables


if __name__ == "__main__":
    # ulx_systems = process_startack.main()
    # ulxlc_runner.main()
    
    # N_lim_dict = create_normalisation_lookup.main()
    # curve_classifications = curve_classifier.main()
    
    # create_erass_curves.main()
    # parent_population = create_erass_parent_population.main()
    
    # erass_simulation.main()
    # process_erass_curves.main()
    print('Hello')


"""

ulx_pop
-------

void        |   
Curve       |   
System      |   
Population  |   populations.py
Table/df    |

"""



# Legacy
# alive_dead_calc.main()
# bh_vs_classification.main()
# eROSITA_lightcurve_sampler
# eROSITA_results_processor