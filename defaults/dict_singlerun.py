# --- DICTIONARIES FOR DEFAULT VALUES
#     Tyler Kukla -- Apr, 2024


# ----------------------------------------------------
# 
# --- SINGLE RUN DICTIONARIES 
# 
# ----------------------------------------------------
singlerun_default = {
    'include_N': True,             # include nitrogen dynamics in solution
    'include_Al': False,           # include aluminum dynamics
    'alphase': 'amal',   # 'gb'    # set aluminum phase
    'use_CaCl2': False,            # relevant for lab run
    'include_DIC': True,           
    'duration': 2,                 # [yr] duration of simulation (or single targetpH or minPH iteration for multi-year)
    'dustsp': 'gbas',              # added dust species 
    'dustsp_2nd': 'amnt',          # added dust species for secondary dust
    'imix': 3,                     # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    'ztot_lab': 0.05,
    'ttot_lab': 100, 
    'water_frac': 1.,              # water fraction for lab simulation
    'spinrun': 'spinup_1_spintuneupManual',         # name of spinup run
    'newrun_id': 'testname',        # name of new run
    'modeldir': '/home/tykukla/SCEPTER/', # directory of model scripts
    'outdir': '/home/tykukla/SCEPTER/scepter_output/', # directory of output files
    'spindir': 's3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/', # location of spinup run
    'climatedir': '/home/tykukla/aglime-swap-cdr/scepter/clim-inputs/',   # climate input main directory
    'climatefiles': 'default',     # climate input subdirectory (contains the climate `.in` files)
    'dustrate': 100,              # [g m-2 yr-1] dust application flux ; divide by 100 to get ton / ha / yr; not used for spintuneups and an initial guess for target ph runs
    'dustrate_2nd': 0,            # [g m-2 yr-1] dust application flux ; multiply by ~11 to get lbs / acre / yr; divide by 100 to get ton / ha / yr (only used for reApp scripts)
    'taudust': 0.05,               # [yr] duration of dust application in year
    'duststart': 0.25,             # [yr] (only used in seasonal runs with v1.0.2 or greater) time of year when dust application starts
    'dustrad': 150,                # [micron] radius of dust particles (gets converted to meters in python script that runs scepter)
    'singlerun_seasonality': True, # [True, False, "spinvalue"] whether to impose seasonality in the singlerun script
    'include_roughness_sa': True,  # [True, False, "spinvalue"] whether to include roughness in the mineral surface area calculation (for switches.in)
    'sa_rule1': "spinvalue",       # [True, False, "spinvalue"] SA decreases as porosity increases
    'sa_rule2': "spinvalue",       # [True, False, "spinvalue"] SA increases as porosity increases
    'include_psd_bulk': True,      # [True, False, "spinvalue"] whether to compute bulk particle size diameters (for switches.in)
    'include_psd_full': True,      # [True, False, "spinvalue"] whether to compute full particle size diameters (for switches.in)
    'poro_iter_field': False,      # [True, False "spinvalue"] porosity iteration
    'poro_evol': False,            # [True, False "spinvalue"] porosity evolves with solid phase dissolution
    'cec_adsorption_on': False,    # [True, False "spinvalue"] whether to enable cec adsorption
    'dep_sample': 0.15,            # [cm] depth of sample for comparing pH to target -- not relevant for single run or initial tuneup, but relevant for tuning rock application
    # --- particle size distribution
    'psdrain_datfile': "psdrain_100um.in",  # [] name of psdrain.in file in SCEPTER/data/ to use for the given run
    'use_psdrain_datfile': False,  # [True, False] If true, the psdrain_datfile is copied over and used, otherwise a single peak distribution is constructed from the gaussian parameters below
    'psdrain_meanRad': 5e-6,       # [m] mean radius of the particle size distribution (only used if `include_psd_*` is True)
    'psdrain_log10_sd': 0.2,       # [m] standard deviation of the particle size distribution in log10 (only used if `include_psd_*` is True)
    'psdrain_wt': 1.,              # [ ] weight for the PSD distribution
    
    # --- optional, update CEC from spinup
    'cec_update_from_spinup': False,   # [True, False] whether to update CEC and alpha vars relative to the spinup value (False means no change to cec.in is made)
    'cec': 21.1,                   # [cmol kg-1] [only used if cec_update_from_spinup == True] cation exchange capacity
    'alpha': 2.,                   # [only used if cec_update_from_spinup == True] pH dependence of CEC coefficients (see Appelo, 1994) -- value from Yoshi (pers. comm.) (usually defined in spinup_*.sh)

    # --- postprocess
    "postproc_prof_list": ["all"],  # list of postproc files to convert to .nc. Options are: ["adsorbed_percCEC", "adsorbed_ppm", "adsorbed", "aqueous_total", 
                                    #                                                         "aqueous", "bulksoil", "exchange_total", "gas", "rate", "soil_ph", "solid_sp_saturation", 
                                    #                                                         "solid_volumePercent", "solid_weightPercent", "solid", "specific_surface_area", "surface_area"]

    # --- OPTIONAL (accepts spinup value if None)
    "qrun": None,                  # [m yr-1] water infiltration flux, use this to override spinup only
    "mat": None,                   # [degC] mean annual temperature, use this to override spinup only
    "dust_mixdep": None,           # [m] depth of dust mixing into soil column
    "soilmoisture_surf": None,     # [] soil moisture saturation at the surface of the profile

    # --- identify scepter version from executable name
    "v102_exelist": ["scepter_richards"],   # list of executables for v1.0.1

    # --- spinup check
    'spinup_on': False,            # [bool] if it's a "spinup" initiated from another spinup, then our file naming convention will be different

    # --- compute specific
    'aws_save': "None",              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/",  # where to save at AWS (only used if 'aws_save'=True)

    # --- which scepter executable to use
    'scepter_exec_name': 'scepter'  # ['scepter', 'scepter_rateA', ...]
}


# ----------------------------------------------------
# 
# --- SPINUP DICTIONARIES 
# 
# ----------------------------------------------------
spinup_default = {
    'include_N': True,             # include nitrogen dynamics in solution
    'include_Al': False,           # include aluminum dynamics
    'alphase': 'amal',   # 'gb'    # set aluminum phase
    'use_CaCl2': False,            # relevant for lab run
    'include_DIC': True,           
    'use_local_storage': False,    # only relevant for georgia tech hpc
    'cec': 21.1,                   # [cmol kg-1] cation exchange capacity
    'modeldir': '/home/tykukla/SCEPTER/', # directory of model scripts
    'outdir': '/home/tykukla/SCEPTER/scepter_output/', # directory of output files
    'include_roughness_sa': True,  # [True, False] whether to include roughness in the mineral surface area calculation (for switches.in)
    'cec_adsorption_on': False,    # [True, False] whether to enable cec adsorption
    'dep_sample': 0.15,            # [cm] depth of sample for comparing pH to target -- not relevant for single run or initial tuneup, but relevant for tuning rock application

    # --- tunespin specific
    'activity_on': False,          # [True, False] whether to turn on thermodynamic activity coefficients (in switches.in file) 
    'make_initial_guess': False,   # [True, False] whether to use existing tuned vars to guess the correct value (e.g., if a nearby site is already spun up). Assigning a `runname_guess` will set to True
    'initguess': 'none',           # [] directory name to pull initial parameter guess (value other than 'none' sets 'make_initial_guess' to True)
    'stop_unsuccessful': True,     # [True, False] whether to stop looking after run timeout threshold is passed
    'add_secondary': True,         # [True, False] whether to add secondary precipitates to list of solids to track (defined in array below)
    'sld_track': ["cc", "ka",      # list of minerals whose secondary precipitation to track and output (appended to the sld_list in python script)
                  "gb", "ct", "cabd", "ill", "gps", "mgbd"],
    'liming': False,               # [True, False] if liming, 10 units of ca added (over-writes initial ca value below) 
    'limesp': 'cc',                # [] lime species to apply (CaCO3 or CaO)
    'water_frac_tunespin': 2.5,    # [] water frac for the tunespin (set separately because default files had 2.5 for tuneup, 1 for basalt)
    'iter_max': 3000,              # [] max iterations to find converged solution
    'tph': 7.2,                   # [] target pH for tunespin (usually defined in spinup_*.sh or, for rock app, the input .csv file)
    'tec': 20.9,                   # [%CEC] target exchangeable acidity (acidsat) (usually defined in spinup_*.sh)
    'tsom': 2.05,                  # [wt%] target soil organic matter (usually defined in spinup_*.sh)
    'tsoilco2': -1.804,            # [log10 atm] target soil pco2 (usually defined in spinup_*.sh)
    'poro': 0.447,                 # [] field porosity (usually defined in spinup_*.sh)
    'soilmoisture': 0.282,         # [m3 / m3] soil moisture (usually defined in spinup_*.sh)
    'alpha': 2.,                   # pH dependence of CEC coefficients (see Appelo, 1994) -- value from Yoshi (pers. comm.) (usually defined in spinup_*.sh)
    'ca': 500e-6,                  # [M] initial calcium at surface (500e-6 should be negligible.. 10 would approximate historical liming practice)
    'mat': 8.23,                   # [degC] mean annual temperature (usually defined in spinup_*.sh)
    'spinname': 'default_UPDATEME', #  [] name of spinuprun, should be overwritten by whatever's defined in spinup_*.sh
    'erosion': 0.001013,           # [m/yr] erosion rate (usually defined in spinup_*.sh)
    'qrun': 0.351,                 # [m/yr] mean annual runoff rate (usually defined in spinup_*.sh)
    'nitrif': 1.005952,            # [gN/m2/yr] NO3 production rate via nitrification (usually defined in spinup_*.sh) (24.6 ~ 220 lbs/acre/year)
    'phnorm_pw': False,            # metric for target pH (true = porewater pH; false = soil pH) -- not relevant for single run
    'spinup_parentrock_file': 'parentlist1.json',  # [None or *.json] looks for a parent rock composition in SCEPTER/data/*.json. None assumes fully inert
    'tunespin_timeout': 60*60*48,      # [s] time before we give up on a tunespin run (7200 is 120*60 = 2 hours)

    # --- postprocess
    "postproc_prof_list": ["all"],  # list of postproc files to convert to .nc. Options are: ["adsorbed_percCEC", "adsorbed_ppm", "adsorbed", "aqueous_total", 
                                    #                                                         "aqueous", "bulksoil", "exchange_total", "gas", "rate", "soil_ph", "solid_sp_saturation", 
                                    #                                                         "solid_volumePercent", "solid_weightPercent", "solid", "specific_surface_area", "surface_area"]

    # --- OPTIONAL INPUTS:
    # defaults are stored in the function and used if the value is None
    'ttot_field': None,            # [yr; default=10000] duration of the spinup field run 
    'ztot_field': None,            # [m; default=0.5] depth of soil column for field run
    'nz': None,                    # [n; default=30] grid cells across ztot_field depth
    'zom': None,                   # [m; default=0.25] mixed layer depth for organic matter mixing
    'omrain_field': None,          # [g C/m2/yr; default=900] organic matter rain rate; initial guess for tunespin_*_newton scripts
    'zwater': None,                # [m; default=10000] depth of water table
    'w_scheme_field': None,        # [default=1] erosion scheme; 0-- cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible(cnst porosity prof), if not defined 0 is taken
    'mix_scheme_field': None,      # [default=1] mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken
    'poro_iter_field': None,       # [default='false'] porosity iteration
    'poro_evol': None,             # [default='false'] porosity evolves with solid phase dissolution
    'sldmin_lim': None,            # [default='false'] limiting mineral lowest concentration (for numerical solver reasons)
    'psd_bulk_field': None,        # [default='false'] enabling PSD tracking
    'psd_full_field': None,        # [default='false'] enabling PSD tracking for individual solid species
    'sa_evol_1': None,             # [default='true'] surface area decreases with increased porosity
    'sa_evol_2': None,             # [default='false'] surface area decreases with increased porosity  
    'display': None,               # [default='true'] display model outputs as it runs
    'disp_lim': None,              # [default='true'] read out only limited model results
    'close_aq_field': None,        # [default='false'] force closed system conditions for aqueous phases
    'season': None,                # [default='false'] allow seasonal clim variability

    # --- identify scepter version from executable name
    "v102_exelist": ["scepter_richards"],   # list of executables for v1.0.1

    # --- compute specific
    'aws_save': "copy",              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/",  # where to save at AWS (only used if 'aws_save'=True)

    # --- which scepter executable to use
    'scepter_exec_name': 'scepter'  # ['scepter', 'scepter_rateA', ...]
}



# ----------------------------------------------------
# 
# --- RE-APPLICATION DICTIONARIES 
# 
# ----------------------------------------------------
# Multi-year type 1 -- specifiy the dust timeseries
specifydust_default = { 
    'include_N': True,             # include nitrogen dynamics in solution
    'include_Al': False,           # include aluminum dynamics
    'alphase': 'amal',   # 'gb'    # set aluminum phase
    'use_CaCl2': False,            # relevant for lab run
    'include_DIC': True,           
    'use_local_storage': False,    # only relevant for georgia tech hpc
    # 'cec': 21.1,                 # [cmol kg-1] cation exchange capacity
    'duration': 1,                 # [yr] duration of simulation (or iteration for multiyear; gets over-written by the dust timeseries file)
    'dustsp': 'gbas',              # added dust species (will get overwritten by dust timeseries file if defined there) 
    'dustsp_2nd': 'amnt',          # added dust species for secondary dust (will get overwritten by dust timeseries file if defined there) 
    'dustrate': 1000,              # [g m-2 yr-1] dust application flux ; divide by 100 to get ton / ha / yr; (will get overwritten by dust timeseries file if defined there) 
    'dustrate_2nd': 80,            # [g m-2 yr-1] dust application flux ; multiply by ~11 to get lbs / acre / yr; divide by 100 to get ton / ha / yr (will get overwritten by dust timeseries file if defined there) 
    'imix': 3,                     # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    # 'catlist': ['ca','mg','k','na'], # list of tracked cations
    'spinrun': 'spinup_1_spintuneupManual',         # name of spinup run
    'newrun_id': 'testnameMeanAnn_gbas_fert',        # name of new run
    'modeldir': '/home/tykukla/SCEPTER/', # directory of model scripts
    'outdir': '/home/tykukla/SCEPTER/scepter_output/', # directory of output files
    'spindir': 's3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/', # location of spinup run
    'climatedir': '/home/tykukla/aglime-swap-cdr/scepter/clim-inputs/',   # climate input main directory
    'climatefiles': 'default',     # climate input subdirectory (contains the climate `.in` files)
    'taudust': 0.05,               # [yr] duration of dust application in year
    'duststart': 0.25,             # [yr] (only used in seasonal runs with v1.0.2 or greater) time of year when dust application starts
    'dustrad': 150,                # [micron] radius of dust particles (gets converted to meters in python script that runs scepter)
    # 'add_secondary': False,         # [True, False] whether to add secondary precipitates to list of solids to track (defined in array below)
    # 'sld_track': ["cc", "ka",      # list of minerals whose secondary precipitation to track and output (appended to the sld_list in python script)
    #               "gb", "ct", "cabd", "ill", "gps", "mgbd"],
    'singlerun_seasonality': False, # [True, False] whether to impose seasonality in the singlerun script
    'include_roughness_sa': False,  # [True, False] whether to include roughness in the mineral surface area calculation (for switches.in)
    'include_psd_bulk': False,      # [True, False] whether to compute bulk particle size diameters (for switches.in)
    'include_psd_full': False,      # [True, False] whether to compute full particle size diameters (for switches.in)
    'cec_adsorption_on': False,    # [True, False] whether to enable cec adsorption
    'dep_sample': 0.15,            # [cm] depth of sample for comparing pH to target -- not relevant for single run or initial tuneup, but relevant for tuning rock application
    # --- particle size distribution
    'psdrain_datfile': "psdrain_100um.in",  # [] name of psdrain.in file in SCEPTER/data/ to use for the given run
    'use_psdrain_datfile': False,  # [True, False] If true, the psdrain_datfile is copied over and used, otherwise a single peak distribution is constructed from the gaussian parameters below
    'psdrain_meanRad': 5e-6,       # [m] mean radius of the particle size distribution (only used if `include_psd_*` is True)
    'psdrain_log10_sd': 0.2,       # [m] standard deviation of the particle size distribution in log10 (only used if `include_psd_*` is True)
    'psdrain_wt': 1.,              # [ ] weight for the PSD distribution

    # --- multi-year specific
    "dust_ts_dir": "/home/tykukla/ew-workflows/inputs/scepter/dust",
    "dust_ts_fn": "gbas_15yr_1app_no2nd_001.csv",
    'max_time': 6,                 # [yr] total amount of time to simulate (must be > `duration`, which refers to individual targetpH run)
    'clim_files': ["T_temp.in", "q_temp.in", "Wet_temp.in"],  # names of climate files (only used in multi-year because we have to update each)

    # --- lab vars
    'ztot_lab': 0.05,
    'ttot_lab': 100, 
    'water_frac': 1.,              # water fraction for lab simulation

    # --- postprocess
    "postproc_prof_list": ["all"],  # list of postproc files to convert to .nc. Options are: ["adsorbed_percCEC", "adsorbed_ppm", "adsorbed", "aqueous_total", 
                                    #                                                         "aqueous", "bulksoil", "exchange_total", "gas", "rate", "soil_ph", "solid_sp_saturation", 
                                    #                                                         "solid_volumePercent", "solid_weightPercent", "solid", "specific_surface_area", "surface_area"]

    # --- identify scepter version from executable name
    "v102_exelist": ["scepter_richards"],   # list of executables for v1.0.1
    
    # --- OPTIONAL (accepts spinup value if None)
    "qrun": None,                  # [m yr-1] water infiltration flux, use this to override spinup only
    "mat": None,                   # [degC] mean annual temperature, use this to override spinup only
    
    # --- compute specific
    'aws_save': None,              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/",  # where to save at AWS (only used if 'aws_save'=True)

    # --- which scepter executable to use
    'scepter_exec_name': 'scepter'  # ['scepter', 'scepter_rateA', ...]
}

