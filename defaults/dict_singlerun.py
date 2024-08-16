# --- DICTIONARIES FOR DEFAULT SINGLERUN VALUES
#     Tyler Kukla -- Apr, 2024

default_dictionary = {
    'include_N': True,             # include nitrogen dynamics in solution
    'include_Al': False,           # include aluminum dynamics
    'alphase': 'amal',   # 'gb'    # set aluminum phase
    'use_CaCl2': False,            # relevant for lab run
    'include_DIC': True,           
    'use_local_storage': False,    # only relevant for georgia tech hpc
    'cec': 21.1,                   # [cmol kg-1] cation exchange capacity
    'duration': 2,                 # [yr] duration of simulation (or single targetpH or minPH iteration for multi-year)
    'dustsp': 'gbas',              # added dust species 
    'dustsp_2nd': 'amnt',          # added dust species for secondary dust
    'imix': 3,                     # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    'ztot_lab': 0.05,
    'ttot_lab': 100, 
    'water_frac': 1.,              # water fraction for lab simulation
    'catlist': ['ca','mg','k','na'], # list of tracked cations
    'spinrun': 'spinup_1_spintuneupManual',         # name of spinup run
    'newrun_id': 'testname',        # name of new run
    'modeldir': '/home/tykukla/SCEPTER/', # directory of model scripts
    'outdir': '/home/tykukla/SCEPTER/scepter_output/', # directory of output files
    'climatedir': '/home/tykukla/aglime-swap-cdr/scepter/clim-inputs/',   # climate input main directory
    'climatefiles': 'default',     # climate input subdirectory (contains the climate `.in` files)
    'dustrate': 100,              # [g m-2 yr-1] dust application flux ; divide by 100 to get ton / ha / yr; not used for spintuneups and an initial guess for target ph runs
    'dustrate_2nd': 0,            # [g m-2 yr-1] dust application flux ; multiply by ~11 to get lbs / acre / yr; divide by 100 to get ton / ha / yr (only used for reApp scripts)
    'taudust': 0.05,               # [yr] duration of dust application in year
    'dustrad': 150,                # [micron] radius of dust particles (gets converted to meters in python script that runs scepter)
    'add_secondary': False,         # [True, False] whether to add secondary precipitates to list of solids to track (defined in array below)
    'sld_track': ["cc", "ka",      # list of minerals whose secondary precipitation to track and output (appended to the sld_list in python script)
                  "gb", "ct", "cabd", "ill", "gps", "mgbd"],
    'singlerun_seasonality': True, # [True, False] whether to impose seasonality in the singlerun script
    'include_roughness_sa': True,  # [True, False] whether to include roughness in the mineral surface area calculation (for switches.in)
    'include_psd_bulk': True,      # [True, False] whether to compute bulk particle size diameters (for switches.in)
    'include_psd_full': True,      # [True, False] whether to compute full particle size diameters (for switches.in)
    
    # --- tunespin specific
    'activity_on': False,          # [True, False] whether to turn on thermodynamic activity coefficients (in switches.in file) 
    'cec_adsorption_on': False,    # [True, False] whether to enable cec adsorption
    'make_initial_guess': False,   # [True, False] whether to use existing tuned vars to guess the correct value (e.g., if a nearby site is already spun up). Assigning a `runname_guess` will set to True
    'stop_unsuccessful': True,     # [True, False] whether to stop looking after run timeout threshold is passed
    'liming': False,               # [True, False] if liming, 10 units of ca added (over-writes initial ca value below) 
    'limesp': 'cc',                # [] lime species to apply (CaCO3 or CaO)
    'limewt': 100.089,             # [] lime wt% of total dust
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
    'initguess': 'none',           # [] directory name to pull initial parameter guess (value other than 'none' sets 'make_initial_guess' to True)
    'spinname': 'default_UPDATEME', #  [] name of spinuprun, should be overwritten by whatever's defined in spinup_*.sh
    'erosion': 0.001013,           # [m/yr] erosion rate (usually defined in spinup_*.sh)
    'qrun': 0.351,                 # [m/yr] mean annual runoff rate (usually defined in spinup_*.sh)
    'nitrif': 1.005952,            # [gN/m2/yr] NO3 production rate via nitrification (usually defined in spinup_*.sh) (24.6 ~ 220 lbs/acre/year)
    'dep_sample': 0.15,            # [cm] depth of sample for comparing pH to target -- not relevant for single run or initial tuneup, but relevant for tuning rock application
    'phnorm_pw': False,            # metric for target pH (true = porewater pH; false = soil pH) -- not relevant for single run
    'maxiter': 50,                 # max iterations for reaching target solution -- not relevant for single run
    'spinup_parentrock_file': 'parentlist1.json',  # [None or *.json] looks for a parent rock composition in SCEPTER/data/*.json. None assumes fully inert
    'tunespin_timeout': 60*60*48,      # [s] time before we give up on a tunespin run (7200 is 120*60 = 2 hours)
    
    # --- multi-year specific
    'max_time': 11,                 # [yr] total amount of time to simulate (must be > `duration`, which refers to individual targetpH run)
    'next_dustrate': 4,            # [g m-2 yr-1] dust rate at the start of the second iteration
    'clim_files': ["T_temp.in", "q_temp.in", "Wet_temp.in"],  # names of climate files (only used in multi-year because we have to update each)

    # --- compute specific
    'aws_save': "None",              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/",  # where to save at AWS (only used if 'aws_save'=True)
    }



s311_dictionary = {
    'alphase': 'amal',   # 'gb'    # set aluminum phase
    'use_CaCl2': False,            # relevant for lab run
    'include_DIC': True,           
    'use_local_storage': False,    # only relevant for georgia tech hpc
    'cec': 21.1,                   # [cmol kg-1] cation exchange capacity
    'duration': 5,                 # [yr] duration of simulation (or single targetpH iteration for multi-year)
    'dustsp': 'gbas',              # added dust species 
    'imix': 3,                     # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    'ztot_lab': 0.05,
    'ttot_lab': 100, 
    'water_frac': 1.,              # water fraction for lab simulation
    'catlist': ['ca','mg','k','na'], # list of tracked cations
    'spinrun': 'site_311',         # name of spinup run
    'newrun_id': 'oldname',        # name of new run
    'modeldir': '/home/tykukla/SCEPTER/', # directory of model scripts
    'outdir': '/home/tykukla/SCEPTER/scepter_output/', # directory of output files
    'climatedir': '/home/tykukla/aglime-swap-cdr/scepter/clim-inputs/',   # climate input main directory
    'climatefiles': 'default',     # climate input subdirectory (contains the climate `.in` files)
    'dustrate': 5000,              # [g m-2 yr-1] dust application flux ; divide by 100 to get ton / ha / yr; not used for spintuneups and an initial guess for target ph runs
    'taudust': 0.05,               # [yr] duration of dust application in year
    'dustrad': 150,                # [micron] radius of dust particles (gets converted to meters in python script that runs scepter)
    
    # --- tunespin specific
    'activity_on': False,          # [True, False] whether to turn on thermodynamic activity coefficients (in switches.in file) 
    'cec_adsorption_on': True,    # [True, False] whether to enable cec adsorption
    'make_initial_guess': False,   # [True, False] whether to use existing tuned vars to guess the correct value (e.g., if a nearby site is already spun up). Assigning a `runname_guess` will set to True
    'stop_unsuccessful': True,     # [True, False] whether to stop looking after run timeout threshold is passed
    'liming': False,               # [True, False] if liming, 10 units of ca added (over-writes initial ca value below) 
    'limesp': 'cc',                # [] lime species to apply (CaCO3 or CaO)
    'limewt': 100.089,             # [] lime wt% of total dust
    'water_frac_tunespin': 2.5,    # [] water frac for the tunespin (set separately because default files had 2.5 for tuneup, 1 for basalt)
    'iter_max': 3000,              # [] max iterations to find converged solution
    'tph': 6.06,                   # [] target pH for tunespin (usually defined in spinup_*.sh or, for rock app, the input .csv file)
    'tec': 20.9,                   # [%CEC] target exchangeable acidity (acidsat) (usually defined in spinup_*.sh)
    'tsom': 2.05,                  # [wt%] target soil organic matter (usually defined in spinup_*.sh)
    'tsoilco2': -1.804,            # [log10 atm] target soil pco2 (usually defined in spinup_*.sh)
    'poro': 0.5,                   # [] field porosity (usually defined in spinup_*.sh)
    'soilmoisture': 0.282,         # [m3 / m3] soil moisture (usually defined in spinup_*.sh)
    'alpha': 2.,                   # pH dependence of CEC coefficients (see Appelo, 1994) -- value from Yoshi (pers. comm.) (usually defined in spinup_*.sh)
    'ca': 500e-9,                  # [M] initial calcium at surface (500e-6 should be negligible.. 10 would approximate historical liming practice)
    'mat': 8.23,                   # [degC] mean annual temperature (usually defined in spinup_*.sh)
    'initguess': 'none',           # [] directory name to pull initial parameter guess (value other than 'none' sets 'make_initial_guess' to True)
    'spinname': 'spinup_1g2zzzY',   #  [] name of spinuprun, should be overwritten by whatever's defined in spinup_*.sh
    'nitrif': 1.005952,            # [gN/m2/yr] NO3 production rate via nitrification (usually defined in spinup_*.sh) (24.6 ~ 220 lbs/acre/year)
    'dep_sample': 0.15,            # [cm] depth of sample for comparing pH to target -- not relevant for single run or initial tuneup, but relevant for tuning rock application
    'phnorm_pw': False,            # metric for target pH (true = porewater pH; false = soil pH) -- not relevant for single run
    'maxiter': 50,                 # max iterations for reaching target solution -- not relevant for single run
    'spinup_parentrock_file': 'parentlist7.json',  # [None or *.json] looks for a parent rock composition in SCEPTER/data/*.json. None assumes fully inert
    'tunespin_timeout': 60*60*48,      # [s] time before we give up on a tunespin run (7200 is 120*60 = 2 hours)

    # --- manual tune spinup specific
    'zom': 0.25,                   # [m] mixed layer depth for organic matter mixing
    'omrain_field': 15.0,           # [g C/m2/yr] organic matter rain rate
    'g2_turnoverTime': 4000,        # [yr] turnover time for g2 class of organic matter (applies to kinsp.in file)
    'zwater': 1000,                # [m] depth of water table
    'erosion': 0.00008692,         # [m/yr] erosion rate (usually defined in spinup_*.sh)
    'qrun': 0.422,                 # [m/yr] mean annual runoff rate (usually defined in spinup_*.sh)
    'w_scheme_field': 1,           # erosion scheme; 0-- cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible(cnst porosity prof), if not defined 0 is taken
    'mix_scheme_field': 1,         # mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken
    'poro_iter_field': 'false',    # porosity iteration
    'sldmin_lim': 'false',         # limiting mineral lowest concentration
    'psd_bulk_field': 'false',     # enabling PSD tracking
    'psd_full_field': 'false',     # enabling PSD tracking for individual solid species
    'include_N': True,             # include nitrogen dynamics in solution
    'include_Al': False,           # include aluminum dynamics
    'add_secondary': False,         # [True, False] whether to add secondary precipitates to list of solids to track (defined in array below)
    'sld_track': ["cc", "ka",      # list of minerals whose secondary precipitation to track and output (appended to the sld_list in python script)
                  "gb", "ct", "cabd", "ill", "gps", "mgbd"],
    'singlerun_seasonality': False,# [True, False] whether to impose seasonality in the singlerun script
    'include_roughness_sa': False,  # [True, False] whether to include roughness in the mineral surface area calculation (for switches.in)

    # --- multi-year specific
    'max_time': 5,                 # [yr] total amount of time to simulate (must be > `duration`, which refers to individual targetpH run)
    'next_dustrate': 5,            # [g m-2 yr-1] dust rate at the start of the second iteration
    'clim_files': ["T_temp.in", "q_temp.in", "Wet_temp.in"],  # names of climate files (only used in multi-year because we have to update each)

    # --- compute specific
    'aws_save': "copy",              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/",  # where to save at AWS (only used if 'aws_save'=True)
    }


reApp_dictionary = {
    'include_N': True,             # include nitrogen dynamics in solution
    'include_Al': False,           # include aluminum dynamics
    'alphase': 'amal',   # 'gb'    # set aluminum phase
    'use_CaCl2': False,            # relevant for lab run
    'include_DIC': True,           
    'use_local_storage': False,    # only relevant for georgia tech hpc
    'cec': 21.1,                   # [cmol kg-1] cation exchange capacity
    'duration': 1,                 # [yr] duration of simulation (or single targetpH or minPH iteration for multi-year)
    'dustsp': 'cc',              # added dust species 
    'dustsp_2nd': 'amnt',            # added dust species for secondary dust
    'imix': 3,                     # mixing style (1=fickian; 2=homogeneous; 3=tilling)
    'ztot_lab': 0.05,
    'ttot_lab': 100, 
    'water_frac': 1.,              # water fraction for lab simulation
    'catlist': ['ca','mg','k','na'], # list of tracked cations
    'spinrun': 'spinup_1_spintuneupManual',         # name of spinup run
    'newrun_id': 'testnameMeanAnn_cc_fert',        # name of new run
    'modeldir': '/home/tykukla/SCEPTER/', # directory of model scripts
    'outdir': '/home/tykukla/SCEPTER/scepter_output/', # directory of output files
    'climatedir': '/home/tykukla/aglime-swap-cdr/scepter/clim-inputs/',   # climate input main directory
    'climatefiles': 'default',     # climate input subdirectory (contains the climate `.in` files)
    'dustrate': 100,               # [g m-2 yr-1] dust application flux ; divide by 100 to get ton / ha / yr; not used for spintuneups and an initial guess for target ph runs
    'dustrate_2nd': 25,            # [g m-2 yr-1] dust application flux ; multiply by ~11 to get lbs / acre / yr; divide by 100 to get ton / ha / yr (only used for reApp scripts)
    'taudust': 0.05,               # [yr] duration of dust application in year
    'dustrad': 150,                # [micron] radius of dust particles (gets converted to meters in python script that runs scepter)
    'add_secondary': False,         # [True, False] whether to add secondary precipitates to list of solids to track (defined in array below)
    'sld_track': ["cc", "ka",      # list of minerals whose secondary precipitation to track and output (appended to the sld_list in python script)
                  "gb", "ct", "cabd", "ill", "gps", "mgbd"],
    'singlerun_seasonality': False, # [True, False] whether to impose seasonality in the singlerun script
    'include_roughness_sa': True,  # [True, False] whether to include roughness in the mineral surface area calculation (for switches.in)
    
    # --- tunespin specific
    'activity_on': False,          # [True, False] whether to turn on thermodynamic activity coefficients (in switches.in file) 
    'cec_adsorption_on': False,    # [True, False] whether to enable cec adsorption
    'make_initial_guess': False,   # [True, False] whether to use existing tuned vars to guess the correct value (e.g., if a nearby site is already spun up). Assigning a `runname_guess` will set to True
    'stop_unsuccessful': True,     # [True, False] whether to stop looking after run timeout threshold is passed
    'liming': False,               # [True, False] if liming, 10 units of ca added (over-writes initial ca value below) 
    'limesp': 'cc',                # [] lime species to apply (CaCO3 or CaO)
    'limewt': 100.089,             # [] lime wt% of total dust
    'water_frac_tunespin': 2.5,    # [] water frac for the tunespin (set separately because default files had 2.5 for tuneup, 1 for basalt)
    'iter_max': 3000,              # [] max iterations to find converged solution
    'tph': 6.8,                   # [] target pH for tunespin (usually defined in spinup_*.sh or, for rock app, the input .csv file)
    'tec': 20.9,                   # [%CEC] target exchangeable acidity (acidsat) (usually defined in spinup_*.sh)
    'tsom': 2.05,                  # [wt%] target soil organic matter (usually defined in spinup_*.sh)
    'tsoilco2': -1.804,            # [log10 atm] target soil pco2 (usually defined in spinup_*.sh)
    'poro': 0.447,                 # [] field porosity (usually defined in spinup_*.sh)
    'soilmoisture': 0.282,         # [m3 / m3] soil moisture (usually defined in spinup_*.sh)
    'alpha': 2.,                   # pH dependence of CEC coefficients (see Appelo, 1994) -- value from Yoshi (pers. comm.) (usually defined in spinup_*.sh)
    'ca': 500e-6,                  # [M] initial calcium at surface (500e-6 should be negligible.. 10 would approximate historical liming practice)
    'mat': 8.23,                   # [degC] mean annual temperature (usually defined in spinup_*.sh)
    'initguess': 'none',           # [] directory name to pull initial parameter guess (value other than 'none' sets 'make_initial_guess' to True)
    'spinname': 'default_UPDATEME', #  [] name of spinuprun, should be overwritten by whatever's defined in spinup_*.sh
    'erosion': 0.001013,           # [m/yr] erosion rate (usually defined in spinup_*.sh)
    'qrun': 0.351,                 # [m/yr] mean annual runoff rate (usually defined in spinup_*.sh)
    'nitrif': 1.005952,            # [gN/m2/yr] NO3 production rate via nitrification (usually defined in spinup_*.sh) (24.6 ~ 220 lbs/acre/year)
    'dep_sample': 0.15,            # [cm] depth of sample for comparing pH to target -- not relevant for single run or initial tuneup, but relevant for tuning rock application
    'phnorm_pw': False,            # metric for target pH (true = porewater pH; false = soil pH) -- not relevant for single run
    'maxiter': 50,                 # max iterations for reaching target solution -- not relevant for single run
    'spinup_parentrock_file': 'parentlist7.json',  # [None or *.json] looks for a parent rock composition in SCEPTER/data/*.json. None assumes fully inert
    'tunespin_timeout': 60*60*48,      # [s] time before we give up on a tunespin run (7200 is 120*60 = 2 hours)
    
    # --- multi-year specific
    'max_time': 6,                 # [yr] total amount of time to simulate (must be > `duration`, which refers to individual targetpH run)
    'next_dustrate': 4,            # [g m-2 yr-1] dust rate at the start of the second iteration
    'clim_files': ["T_temp.in", "q_temp.in", "Wet_temp.in"],  # names of climate files (only used in multi-year because we have to update each)

    # --- compute specific
    'aws_save': "move",              # ["move", "copy", None] whether to "move" file to aws, just copy it, or nothing at all
    'aws_bucket': "s3://carbonplan-carbon-removal/SCEPTER/scepter_output_scratch/",  # where to save at AWS (only used if 'aws_save'=True)
    }


