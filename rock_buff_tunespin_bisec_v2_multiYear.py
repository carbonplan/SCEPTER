import os
import re
import sys
import time
import numpy as np
import shutil
import get_int_prof
import make_inputs
import random
import defaults.dict_singlerun

# -------------------------------------------------------------
# Function to parse arguments
def parse_arguments(args):
    parsed_args = {}
    i = 1  # Start from index 1 to skip the script name (sys.argv[0])
    while i < len(args):
        if args[i].startswith('--'):
            key = args[i][2:]  # Remove '--' prefix
            if i + 1 < len(args) and not args[i + 1].startswith('--'):
                value = args[i + 1]
                parsed_args[key] = value
                i += 1  # Skip the next item as it's the value for the current key
            else:
                parsed_args[key] = None  # If no value provided, set to None
        i += 1
    return parsed_args

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to set global variables from defaults / system args
def set_vars(default_args, system_args):
    # define pattern for identifying floats in sys.args
    float_pattern = r'^[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?$'
    # dict to save
    save_vars = {}
    for key, value in default_args.items():
        if key in system_args:
            float_test1 = re.match(float_pattern, system_args[key]) is not None # (captures all cases but "2.")
            float_test2 =  system_args[key].replace('.', '', 1).isdigit()  # (misses negatives but gets others including "2.")
            if float_test1 or float_test2:  # check is sys arg is a float
                save_vars[key] = float(system_args[key])
                globals()[key] = float(system_args[key])
            else:
                save_vars[key] = system_args[key]
                globals()[key] = system_args[key]
        else:
            save_vars[key] = value
            globals()[key] = value
    return save_vars
# Function to save the combined dictionary to the run dir
def save_dict_to_text_file(dictionary, filename, delimiter='\t'):
    with open(filename, 'w') as file:
        file.write(f"*** variables set by dictionary and system args\n")
        file.write(f"    (note not all vars are used!!)\n")
        for key, value in dictionary.items():
            file.write(f"{key}{delimiter}{value}\n")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to copy climate variables
def copy_files(src_dir, dst_dir):
    # directory must exist
    if not os.path.exists(dst_dir):
        print("Cannot find " + dst_dir)
    
    # iterate over source files
    for filename in os.listdir(src_dir):
        src_file = os.path.join(src_dir, filename)
        dst_file = os.path.join(dst_dir, filename)
        
        # copy them over (skipping any directories)
        if os.path.isfile(src_file):
            shutil.copy2(src_file, dst_file)
            # print(f"Copied {src_file} to {dst_file}")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to generate timesteps for multi-year
def generate_timesteps(total_duration, timestep):
    timesteparr = []
    current_time = 0
    while current_time <= (total_duration - timestep):
        timesteparr.append(current_time)
        current_time += timestep
    return timesteparr

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to write file marking the current iteration
def write_iter_file_with_marker(values, current_index, filename):
    with open(filename, 'w') as file:
        file.write(f"year \n")
        for i, value in enumerate(values):
            if i == current_index:
                file.write(f"{value} <--- \n")
            else:
                file.write(f"{value}\n")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# function to update climate vars to current iteration
def update_clim(inputfile, outputfile, timezero):
    # open the input and output files
    with open(inputfile, 'r') as f_in, open(outputfile, 'w') as f_out:
        # skip the header line and copy it to new file
        header = next(f_in)
        f_out.write(header)  # write to output
        # read the file line by line
        for line in f_in:
            # split each line into year and climate value
            year, clim = line.strip().split('\t')  # Adjust the delimiter as needed
    
            # convert year to an integer and subtract 5
            year = float(year) - timezero
            
            # check if the adjusted year is greater than or equal to tstep
            if year >= 0:
                formatted_yr = "{:.7f}".format(float(year))
                formatted_clim = "{:.7f}".format(float(clim))
                # write the adjusted year and temperature to the output file
                # print(f"{formatted_yr}\t{formatted_clim}\n")
                f_out.write(f"{formatted_yr}\t{formatted_clim}\n")  # Adjust the delimiter as needed
                
# -------------------------------------------------------------
# -------------------------------------------------------------
# --- set default and system args
sys_args = parse_arguments(sys.argv)   # parse system args
import_dict = sys_args['default_dict'] # set the dictionary to use from system args
def_args = getattr(defaults.dict_singlerun, import_dict)  # get dict attribute

# set global variables
combined_dict = set_vars(def_args, sys_args)  # default unless defined in sys_args

targetpH = tph
tau = duration
added_sp = dustsp
spinid = spinrun
fdust = dustrate
taudust=taudust
dustrad = float(dustrad)/1e6  #  /1e6 converts micron to meters

expid = newrun_id

outdir = outdir
datadir = os.path.join(modeldir, 'data/')


# -------------------------------------------------------------------------------
# LOOP THROUGH TIME STEPS 

# ... get time step array
timestep_dur = duration  # just for clarity... duration is a timestep here
mytsteps = generate_timesteps(max_time, timestep_dur)
counter = 0  # for tracking and file naming (MUST START AT ZERO)
runname_field_old, runname_lab_old = "placeholder1", "placeholder2"   # placeholders for update later

for tstep in mytsteps:
    
    # set the spinup 
    if counter == 0:  # first step uses a set spinup run
        spinup_field    = spinid+'_spintuneup_field'  # this gets over-written with later iterations
        spinup_lab      = spinid+'_spintuneup_lab'    # this gets over-written with later iterations
    else:
        spinup_field = runname_field_old
        spinup_lab = runname_lab_old
        fdust = next_dustrate  # set fdust to some ideally smaller value for the successive iterations
    
    # set runnames
    runname_field   = f"{expid}_startyear-{str(tstep).replace('.','p')}_iter-{int(counter)}_field"
    runname_lab     = f"{expid}_startyear-{str(tstep).replace('.','p')}_iter-{int(counter)}_lab"

    # duplicate directories from spinups
    for (runname,spinup) in [(runname_field,spinup_field),(runname_lab,spinup_lab)]:
        
        src = os.path.join(outdir, spinup)
        dst = os.path.join(outdir, runname)
    
        if not os.path.exists(dst): 
            shutil.copytree(src, dst)
        else:
            shutil.rmtree(dst)
            shutil.copytree(src, dst)

        # save file denoting the iteration
        fn_itermarker = os.path.join(dst, "multiyear-iter.res")
        write_iter_file_with_marker(mytsteps, counter, fn_itermarker)
        # save file denoting the variables used
        combined_dict['dustrate'] = fdust # update dust rate
        fn_dict_save = os.path.join(dst, "vars.res")
        save_dict_to_text_file(combined_dict, fn_dict_save, delimiter='\t')

    
    # duplicate the climate files
    for runname in [runname_field,runname_lab]:
        for thisfile in clim_files:  # loop through all three climate inputs
            src_clim = os.path.join(climatedir, climatefiles, thisfile)
            dst_clim = os.path.join(outdir, runname, thisfile)
            # read from source, update, save to dst
            update_clim(src_clim, dst_clim, tstep)


    
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    # BELOW IS COPIED FROM THE SINGLE YEAR VERSION
    error = 1e4
    tol = 1e-4
    
    exename = 'scepter'
    # exename_src = 'scepter_test'
    exename_src = 'scepter'
    to = ' '
    where = '/'
    
    # ============ common input file modification wrt spinup: for field run
    filename = '/switches.in'
    src = outdir + spinup_field + filename
    dst = outdir + runname_field  + filename
    with open(src, 'r') as file:
        data = file.readlines()
    data[2] = '{:d}\tbio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken\n'.format(imix)
    data[7] = 'true\trestart from a previous run\n'
    data[-3] = 'true\tenabling PSD tracking\n'
    data[-2] = 'true\tenabling PSD tracking for individual solid species\n'
    data[-1] = 'true\tenabling seasonality\n'   # added per yoshi suggestion
    with open(dst, 'w') as file:
        file.writelines(data)
    
    
    if added_sp == 'gbas': dustsrc = os.path.join(modeldir, 'data', 'dust_gbasalt.in')
    if added_sp == 'cc': dustsrc = os.path.join(modeldir, 'data', 'dust_lime.in')
    dustdst = 'dust.in'
    
    os.system('cp ' + dustsrc + to + outdir + runname_field + where + dustdst) 
    
    
    filename = '/slds.in'
    src = outdir + spinup_field + filename
    dst = outdir + runname_field  + filename
    with open(src, 'r') as file:
        data = file.readlines()
    
    data.insert(1, added_sp+'\t\n')
    with open(dst, 'w') as file:
        file.writelines(data)
        
    # ============ adding Fe(II) as tracer and its oxidation =================
    # filename = '/solutes.in'
    # src = outdir + spinup + filename
    # dst = outdir + runname  + filename
    # with open(src, 'r') as file:
        # data = file.readlines()
    
    # data.insert(1, 'fe2\t\n')
    # with open(dst, 'w') as file:
        # file.writelines(data)
        
    # filename = '/gases.in'
    # src = outdir + spinup + filename
    # dst = outdir + runname  + filename
    # with open(src, 'r') as file:
        # data = file.readlines()
    
    # data.insert(1, 'po2\t\n')
    # with open(dst, 'w') as file:
        # file.writelines(data)
        
    # filename = '/extrxns.in'
    # src = outdir + spinup + filename
    # dst = outdir + runname  + filename
    # with open(src, 'r') as file:
        # data = file.readlines()
    
    # data.insert(1, 'fe2o2\t\n')
    # with open(dst, 'w') as file:
        # file.writelines(data)
        
    # ============ common input file modification wrt spinup: for lab run ============
    
    filename = '/slds.in'
    src = outdir + spinup_lab + filename
    dst = outdir + runname_lab  + filename
    with open(src, 'r') as file:
        data = file.readlines()
    
    data.insert(1, added_sp+'\t\n')
    with open(dst, 'w') as file:
        file.writelines(data)
        
    # filename = '/solutes.in'
    # src = outdir + spinup_lab + filename
    # dst = outdir + runname_lab  + filename
    # with open(src, 'r') as file:
        # data = file.readlines()
    
    # data.insert(1, 'k\t\nna\t\nmg\t\n')
    # with open(dst, 'w') as file:
        # file.writelines(data)
        
    filename = '/kinspc.in'
    src = outdir + spinup_lab + filename
    dst = outdir + runname_lab  + filename
    with open(src, 'r') as file:
        data = file.readlines()
        
    data.insert(1, added_sp+'\t0\n')
    with open(dst, 'w') as file:
        file.writelines(data)
        
    filename = '2ndslds.in'
    srcfile = os.path.join(datadir, '2ndslds_def.in')
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,srcfile = srcfile
        )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,srcfile = srcfile
        )
        
    # compile 
    # os.system('make')
    # os.system('make --file=makefile_test')
    # for runname in [runname_field,runname_lab]:
        # if not os.path.exists( outdir + runname) : 
            # os.system('mkdir -p ' + outdir + runname)
            # os.system('cp ' + exename_src + to + outdir + runname + where + exename)
    
    maxiter = 50

    lastrun = False  # [TK added] turns to true if special condition requires ending run before error < tol
    
    res_list = []
    cnt = 0
    
    # get ph of spin-up
    phint = get_int_prof.get_ph_int_site(outdir,spinup_lab,dep_sample)
    phint_field = get_int_prof.get_ph_int_site(outdir,spinup_field,dep_sample)
    ymx = phint - targetpH
    if phnorm_pw: ymx = phint_field - targetpH
    res_list.append([cnt, phint_field,phint, targetpH, 0, abs( ymx/targetpH ) ])
    
    cnt += 1
    
    while (error > tol):
        
        
        ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ## --- setup for field run --- ##
        
        filename = '/frame.in'
        src = outdir + spinup_field + filename
        dst = outdir + runname_field  + filename
    
        with open(src, 'r') as file:
            data = file.readlines()
        data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(tau)
        data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust)
        data[7]     = '{:.8f}\tduration of dust application [yr]\n'.format(taudust)
        data[16]    = '{:.8f}\tradius of particles [m]\n'.format(dustrad)   # [tykukla added]
        data[18]    = '{}\n'.format(spinup_field)
        data[20]    = '{}\n'.format(runname_field)
        with open(dst, 'w') as file:
            file.writelines(data)
            
            
        ## --- run field run --- ##
        
        print(outdir+runname_field+'/scepter')
        os.system(outdir+runname_field+'/scepter')
    
        ## --- getting data from field run --- ##
        
        phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
        dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
        sps = ['g2','inrt',added_sp]
        if include_Al:  sps.append( alphase )
        sldwt_list = []
        for sp in sps:
            sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
            sldwt_list.append(sldwt)
        
        aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
        print(aqsps,btmconcs,dep)
        
        dic,dep = get_int_prof.get_ave_DIC_save(outdir,runname_field,dep_sample) # returning DIC in mol/ solid m3
        
        ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ## --- setup for lab run --- ##
        
        poro_lab = water_frac/(1./dense_lab+water_frac)
        
        oxide_ctnm_list = ['ca','mg','na','k'] 
        oxide_oxnm_list = ['cao','mgo','na2o','k2o'] 
        oxide_stch_list = [1,1,2,2] 
        oxide_mass_list = [56.1 ,40.3, 62, 94.2]
        
        if include_N:
            oxide_ctnm_list.append( 'no3' )
            oxide_oxnm_list.append( 'amnt' )
            oxide_stch_list.append( 2 )
            oxide_mass_list.append( 80 )
        if include_Al:
            oxide_ctnm_list.append( 'al' )
            oxide_oxnm_list.append( 'al2o3' )
            oxide_stch_list.append( 2 )
            oxide_mass_list.append( 1.02E+02 )
        
        fdust_list = []
        fdust_nm_list = []
    
        for sp in aqsps:
            isp = aqsps.index(sp)
            iox = oxide_ctnm_list.index(sp)
            conc = btmconcs[isp]
            fdust_tmp = ztot_lab*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
            fdust_list.append(fdust_tmp)
            fdust_nm_list.append(oxide_oxnm_list[iox])
            
            
        fdust_lab = fdust_list[aqsps.index('ca')]
        
        if use_CaCl2:
            cacl2_conc = 0.01
            cacl2_wt2  = 110.98
            fdust_cacl2 = ztot_lab*poro_lab*cacl2_conc*1e3*cacl2_wt2
        
            if fdust_cacl2 > 0:
                fdust_list.append(fdust_cacl2)
                fdust_nm_list.append('cacl2')
                
                
        if include_DIC:                                     # (added 3.23.2023)
            fdust_dic = ztot_lab*(1-poro_lab)*dic* 30.      # (added 3.22.2023)
            fdust_list.append(fdust_dic)                    # (added 3.22.2023)
            fdust_nm_list.append('g1')                      # (added 3.22.2023)
    
        fdust_list = [fdust_tmp/fdust_lab  for fdust_tmp in fdust_list  ]  
        
        filename = '/frame.in'
        src = outdir + spinup_lab + filename
        dst = outdir + runname_lab  + filename
    
        with open(src, 'r') as file:
            data = file.readlines()
        data[1]     = '{:.8f}\ttotal depth of weathering profile [m]\n'.format(ztot_lab)
        data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(ttot_lab)
        data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust_lab)
        data[10]     = '{:.8f}\tinitial porosity\n'.format(poro_lab)
        data[16]    = '{:.8f}\tradius of particles [m]\n'.format(dustrad)   # [tykukla added]
        with open(dst, 'w') as file:
            file.writelines(data)
            
            
        filename = 'dust.in'
        sld_varlist = [ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ] 
        make_inputs.get_input_sld_properties(
            outdir=outdir
            ,runname=runname_lab
            ,filename = filename
            ,sld_varlist=sld_varlist
            )
            
        pr_list_lab = [(sp,sldwt_list[sps.index(sp)]/100.) for sp in sps]
        atm_list_lab = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
        if include_DIC: atm_list_lab = [('pco2',1e-20),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)] # (added 3.22.2023)
        make_inputs.get_input_tracer_bounds(
            outdir=outdir
            ,runname=runname_lab
            ,pr_list = pr_list_lab
            # ,rain_list=rain_list
            ,atm_list=atm_list_lab
            )
        
        # break
        
        ## --- run lab run --- ##
        
        print(outdir+runname_lab+'/scepter')
        os.system(outdir+runname_lab+'/scepter')
        
        ## --- getting data from lab run --- ##
        
        phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
        print(phint_field,phint)
        
        ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ## --- Newton + bisection iteration --- ## 
        
        time.sleep(5)
        
    
        ymx = phint - targetpH
        if phnorm_pw: ymx = phint_field - targetpH
        ymx = -ymx
            
    
        emx = ymx/targetpH
        error = np.max(np.abs(emx))
    
        print(fdust)
        
        print(' ')
        print(' ')
        print(' ')
        print(' ')
        print("******* {:d}'s iteration, error = {:.8f}".format(cnt,error))
        print(' ')
        print(' ')
        print(' ')
        print(' ')
        
        res_list.append([cnt, phint_field, phint, targetpH, fdust, abs( ymx/targetpH ) ])
        
        cnt += 1

        # --- [TK ADDED] handle special cases
        if lastrun:
            error = 1e-99  # end the while loop
        # ---
        
        
        if np.isnan(ymx):
            print('nan detected in solution')
            error = 1e-99
        else:
            data_tmp = np.array(res_list)
            iph = 1 
            idust = 3 
            if not phnorm_pw: iph = 2
            idust = 4
            
            # --- [TK added: start] 
            data_tmp_skip1 = data_tmp[1:,:]  # skip first row since it inherits the last iter 
            # ... check for two special cases: (1) near-zero application results in pH that's too high; (2) max application is not max pH (non-linear)
            # get indices for min / max pH and dust app
            phmax_idx = np.argmax(data_tmp_skip1[:,iph])      # index for max pH 
            phmin_idx = np.argmin(data_tmp_skip1[:,iph])      # index for min pH 
            dustmax_idx = np.argmax(data_tmp_skip1[:,idust])  # index for max dust
            dustmin_idx = np.argmin(data_tmp_skip1[:,idust])  # index for min dust

            # ... test condition 1 (pH is too high compared to target
            low_dust_thresh = 0.1   # minimum application we allow
            targetpH_plus_buffer = targetpH + (1.1 * tol)  # allow range w/in tolerance
            # if minimum dust flx is below threshold AND associated pH exceeds target pH, we're maxed out
            ph_maxedout_condition =  (np.min(data_tmp_skip1[dustmin_idx,idust]) <= low_dust_thresh) & (data_tmp_skip1[dustmin_idx,iph] > targetpH_plus_buffer)

            # ... test condition 2 (non-linearity)
            non_linear_condition = phmax_idx != dustmax_idx

            # ... handle conditions
            if ph_maxedout_condition:
                print("pH above target despite no application")
                fdust = 0
                lastrun = True
            elif non_linear_condition:
                print("increasing rock application is not increasing pH (and pH below target)")
                fdust = data_tmp_skip1[phmax_idx, idust]  # take the highest pH dust flux and move on...
                lastrun = True
            # --- [TK added: end]
            
            elif np.max(data_tmp[:,iph]) < targetpH:  # if we've never hit the target (and not non-linear) keep increasing fdust
                fdust = np.max(data_tmp[:,idust])*1.5
                print(np.max(data_tmp[:,iph]), targetpH, fdust)
            else:
                # if cnt >2:  fdust_old_old = fdust_old
                # fdust_old = fdust
                imin = np.argmin(1./(data_tmp[:,iph]-targetpH))
                imax = np.argmax(1./(data_tmp[:,iph]-targetpH))
                # fdust = 0.5 * (data_tmp[imin,idust] +  data_tmp[imax,idust])
                # define slope = DpHmax_min/Dfdustmax_min
                # and Dfmax_x = DpHmax_x/slope
                # and then fx = fmax - Dfmax_x
                slope = (data_tmp[imax,iph] - data_tmp[imin,iph])/(data_tmp[imax,idust] - data_tmp[imin,idust])
                fdust = data_tmp[imax,idust] - ( data_tmp[imax,iph] - targetpH )/slope
                
                # detecting if the updated fdust is already tested
                # if it is, it means non-linearlity
                iclose = np.argmin(np.abs(data_tmp[:,idust]-fdust))
                fdust_old =  data_tmp[iclose,idust]
                
                # just avoid repeated solution by randomly picking up from most likely range 
                while abs((fdust- fdust_old)/fdust) < 1e-6:  # the updated fdust is already tested  
                    fdust = random.uniform(data_tmp[imin,idust],data_tmp[imax,idust])
                    iclose = np.argmin(np.abs(data_tmp[:,idust]-fdust))
                    fdust_old =  data_tmp[iclose,idust]
                
                if abs((fdust- fdust_old)/fdust) < 1e-6:  # the updated fdust is already tested  
                    if (data_tmp[iclose,iph] - targetpH)* (data_tmp[imax,iph] - targetpH) >= 0.: 
                        # if the tested value is above target pH, new slope and approximation is calculated switching imax with iclose
                        # Newton-ish new solution
                        # slope = (data_tmp[iclose,iph] - data_tmp[imin,iph])/(data_tmp[iclose,idust] - data_tmp[imin,idust])
                        # fdust = data_tmp[iclose,idust] - ( data_tmp[iclose,iph] - targetpH )/slope
                        
                        # Bisection-ish new solution
                        fdust = 0.5*(data_tmp[iclose,idust] + data_tmp[imin,idust])
                    elif (data_tmp[iclose,iph] - targetpH)* (data_tmp[imin,iph] - targetpH) >= 0.:
                        # if the tested value is below target pH, new slope and approximation is calculated switching imin with iclose
                        # Newton-ish new solution
                        # slope = (data_tmp[imax,iph] - data_tmp[iclose,iph])/(data_tmp[imax,idust] - data_tmp[iclose,idust])
                        # fdust = data_tmp[imax,idust] - ( data_tmp[imax,iph] - targetpH )/slope
                        
                        # Bisection-ish new solution
                        fdust = 0.5*(data_tmp[iclose,idust] + data_tmp[imax,idust])
                    else:
                        # there must be something is wrong
                        print('there must be something wrong')
                        fdust = np.nan
                    # fdust = fdust_old * 0.5
                
                # if cnt >2 and abs((fdust- fdust_old_old)/fdust) < 1e-6:  fdust = fdust_old_old * 0.5
                print(data_tmp[imin,iph],data_tmp[imax,iph], targetpH 
                    ,data_tmp[imin,idust],data_tmp[imax,idust], fdust)
        
        time.sleep(5)
        
        if cnt > maxiter: break
        
        for runname in [runname_field,runname_lab]:
            np.savetxt(outdir + runname + where + 'iteration_tmp.res',np.array(res_list))
        
        
    
    name_list = [
        'iter.'
        ,'porewater_pH[-]'
        ,'soil_pHw[-]'
        ,'target_pH[-]'
        ,'dust[g/m2/yr]'
        ,'error'
        ]
    
    for runname in [runname_field,runname_lab]:
        dst = outdir + runname + where + 'iteration.res'
    
        with open(dst, 'w') as file:
            for item in name_list:
                if name_list.index(item)==len(name_list)-1:
                    file.write('{}\n'.format(item))
                else:
                    file.write('{}\t'.format(item))
            for j in range(len(res_list)):
                item_list = res_list[j]
                for i in range(len(item_list)):
                    item = item_list[i]
                    if i==0:
                        file.write('{:d}\t'.format(item))
                    elif i==len(item_list)-1:
                        file.write('{:.6e}\n'.format(item))
                    else:
                        file.write('{:.6e}\t'.format(item))
                
        print(res_list)
        
        
    if use_local_storage:
        for runname in [runname_field,runname_lab]:
            src = outdir + runname 
            dst = outdir_src + runname 
            
            if not os.path.exists(dst): 
                shutil.copytree(src, dst)
            else:
                shutil.rmtree(dst)
                shutil.copytree(src, dst)


    # UPDATE THE RUNNAME FOR NEXT ITERATION
    counter += 1
    runname_field_old = runname_field
    runname_lab_old = runname_lab