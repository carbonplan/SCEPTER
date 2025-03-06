# %%
# ---------------------------------------------------------
# 
# multiyear run with specified dust flux timeseries (e.g.,
# if you want to apply rock once then run a bunch of years 
# with no application)
# 
# ---------------------------------------------------------
import os
import re
import sys
import time
import numpy as np
import shutil
import get_int_prof
import make_inputs
import random
import subprocess
import defaults.dict_singlerun
import pandas as pd
# import build_composite_multiyear as cfxns


# --- read in helper functions from aglime-swap-cdr
# add aglime-swap-cdr dir to path # [UPDATE FOR YOUR MACHINE]
sys.path.append(os.path.abspath('/home/tykukla/aglime-swap-cdr/scepter/setup'))
# import module
import scepter_helperFxns as shf
import build_composite_multiyear as cfxns
import cflx_proc as cflx
# ---



# %%
# -------------------------------------------------------------
# --- set default and system args
sys_args = shf.parse_arguments(sys.argv)   # parse system args
import_dict = sys_args['default_dict'] # set the dictionary to use from system args
def_args = getattr(defaults.dict_singlerun, import_dict)  # get dict attribute

# set global variables
combined_dict = shf.set_vars(def_args, sys_args)  # default unless defined in sys_args
# (set as kwargs CHANGE WHEN THIS BECOMES A FUNCTION)
kwargs = combined_dict.copy()
# add to globals (CHANGE/REMOVE WHEN THIS BECOMES A FUNCTION)
for key, value in combined_dict.items():
    globals()[key] = value
    
# rename vars
spinid = spinrun
expid = newrun_id
dustrad_init = dustrad  # make a static version since dustrad is overwritten later
datadir = os.path.join(modeldir, 'data/')

# %% 
# --- Read in the dust flux file
# [2] read in file
src_dust = os.path.join(dust_ts_dir, dust_ts_fn)
df_dust = pd.read_csv(src_dust)
# except:
#     ValueError("For now, `rock_buff_dust-ts_multiyear.py` requires a dust file to be defined. We didn't find one.")



# %%
# -------------------------------------------------------------------------------
# LOOP THROUGH TIME STEPS 

# ... get time step array
counter = 0  # for tracking and file naming (MUST START AT ZERO)
runname_field_old, runname_lab_old = "placeholder1", "placeholder2"   # placeholders for update later
# empty lists for moving output to aws
runname_lab_list, runname_field_list = [], []

# --- files to delete if they are copied over
# (diagnostic files from the spinup run)
files_to_delete = ["check_logs.res", "check_results.res", "completed.res"]

# %% 
# Loop through each row in the df
for index, row in df_dust.iterrows():
    # set up the dust values 
    # ********
    # NOTE update this later to work for any arbitrary column header and value 
    # in the .csv file !! 
    # ********
    # [tau] timestep (years)
    if "duration" in row:
        try:
            tau = float(row['duration'])
        except:
            tau = duration
    else:
        tau = duration
    # [dustsp] primary dust
    if "dustsp" in row:
        added_sp = str(row['dustsp'])
    else:
        added_sp = dustsp
    # [dustsp_2nd] secondary dust
    if "dustsp_2nd" in row:
        added_sp2 = str(row['dustsp_2nd'])
    else:
        added_sp2 = dustsp_2nd
    # [dustrate] how much dust to apply 
    if "dustrate" in row:
        try:
            fdust = float(row['dustrate'])
        except:
            fdust = dustrate
    else: 
        fdust = dustrate
    # [dustrate_2nd] how much dust to apply 
    if "dustrate_2nd" in row:
        try:
            fdust2 = float(row['dustrate_2nd'])
        except:
            fdust2 = dustrate_2nd
    else: 
        fdust2 = dustrate_2nd
    # [dustrad] how much dust to apply 
    if "dustrad" in row:
        try:
            dustrad = float(row['dustrad'])/1e6
        except:
            dustrad = float(dustrad_init)/1e6
    else: 
        dustrad = float(dustrad_init)/1e6
    if "yr_start" in row:
        tstep = row["yr_start"]
    else:
        tstep = "NA"
    # ************************************************************************
    
    # set the spinup 
    if counter == 0:  # first step uses a set spinup run
        spinup_field    = spinid+'_field'  # this gets over-written with later iterations
        spinup_lab      = spinid+'_lab'    # this gets over-written with later iterations
    else:
        spinup_field = runname_field_old
        spinup_lab = runname_lab_old
        # fdust = next_dustrate  # set fdust to some ideally smaller value for the successive iterations
    
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
        shf.write_iter_file_with_marker(df_dust['yr_start'].values, counter, fn_itermarker)
        # save file denoting the variables used
        combined_dict['dustrate'] = fdust # update dust rate
        fn_dict_save = os.path.join(dst, "vars.res")
        shf.save_dict_to_text_file(combined_dict, fn_dict_save, delimiter='\t')
        # delete the diagnistic files from the spinup
        for file in files_to_delete:
            if os.path.exists(os.path.join(dst, file)):
                os.remove(os.path.join(dst, file))
    
    if singlerun_seasonality:   # duplicate the climate files
        for runname in [runname_field,runname_lab]:
            for thisfile in clim_files:  # loop through all three climate inputs
                src_clim = os.path.join(climatedir, climatefiles, thisfile)
                dst_clim = os.path.join(outdir, runname, thisfile)
                # read from source, update, save to dst
                shf.update_clim(src_clim, dst_clim, tstep)
    
    # save the dust file
    for runname in [runname_field,runname_lab]:
        src_dust = os.path.join(dust_ts_dir, dust_ts_fn)
        dst_dust = os.path.join(outdir, runname, dust_ts_fn)
        if os.path.isfile(src_dust):
            shutil.copy2(src_dust, dst_dust)


    
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
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
    data[2] = '{:d}\tbio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken\n'.format(int(imix))
    data[7] = 'true\trestart from a previous run\n'
    if include_psd_bulk:
        data[-3] = 'true\tenabling PSD tracking\n'
    else:
        data[-3] = 'false\tenabling PSD tracking\n'
    if include_psd_full:
        data[-2] = 'true\tenabling PSD tracking for individual solid species\n'
    else:
        data[-2] = 'false\tenabling PSD tracking for individual solid species\n'
    if include_roughness_sa == True:
        data[8] = 'true\tinclude roughness in mineral surface area\n'
    else:
        data[8] = 'false\tinclude roughness in mineral surface area\n'
    if singlerun_seasonality==True:
        data[-1] = 'true\tenabling seasonality\n'   # added per yoshi suggestion
    else:
        data[-1] = 'false\tenabling seasonality\n'
    if cec_adsorption_on == True:
        data[11] = "true\tenabling adsorption for cation exchange\n"
    else:
        data[11] = "false\tenabling adsorption for cation exchange\n"

    with open(dst, 'w') as file:
        file.writelines(data)
    

    # --- PRIMARY DUST FILE
    if added_sp == "amnt": dustsrc = os.path.join(modeldir, 'data', 'dust_fert.in')
    if added_sp == 'gbas': dustsrc = os.path.join(modeldir, 'data', 'dust_gbasalt.in')
    if added_sp == 'cc': dustsrc = os.path.join(modeldir, 'data', 'dust_lime.in')
    if added_sp == 'cao': dustsrc = os.path.join(modeldir, 'data', 'dust_cao.in')
    if added_sp == 'dlm': dustsrc = os.path.join(modeldir, 'data', 'dust_dlm.in')
    if added_sp == 'wls': dustsrc = os.path.join(modeldir, 'data', 'dust_wls.in')
    dustdst = 'dust.in'
    
    os.system('cp ' + dustsrc + to + outdir + runname_field + where + dustdst) 

    # --- SECONDARY DUST FILE
    if added_sp2 == "amnt": dustsrc2 = os.path.join(modeldir, 'data', 'dust_fert.in')
    if added_sp2 == 'gbas': dustsrc2 = os.path.join(modeldir, 'data', 'dust_gbasalt.in')
    if added_sp2 == 'cc': dustsrc2 = os.path.join(modeldir, 'data', 'dust_lime.in')
    if added_sp2 == 'cao': dustsrc2 = os.path.join(modeldir, 'data', 'dust_cao.in')
    if added_sp2 == 'dlm': dustsrc2 = os.path.join(modeldir, 'data', 'dust_dlm.in') 
    if added_sp2 == 'wls': dustsrc = os.path.join(modeldir, 'data', 'dust_wls.in')
    dustdst2 = 'dust_2nd.in'
    
    os.system('cp ' + dustsrc2 + to + outdir + runname_field + where + dustdst2) 
    # --- 
    
    filename = '/slds.in'
    src = outdir + spinup_field + filename
    dst = outdir + runname_field  + filename
    with open(src, 'r') as file:
        data = file.readlines()
    if not data[-1].endswith("\n"):
        data[-1] += "\n"   # add to avoid a messy append
    # add dust sp
    data.insert(1, added_sp+'\n')
    if added_sp2 in ['gbas','cc','cao','dlm','amnt', 'wls']: # then add this as well
        data.insert(1, added_sp2+'\n')
    if include_N and added_sp2 != "amnt":
        data.append('amnt'+'\n')
    if include_Al: 
        data.append(alphase+'\n')

    with open(dst, 'w') as file:
        file.writelines(data)
    # remove duplicate minerals
    shf.remove_duplicates(dst)

        
    # ============ adding Fe(II) as tracer and its oxidation =================
    filename = '/solutes.in'
    src = outdir + spinup + filename
    dst = outdir + runname_field  + filename
    with open(src, 'r') as file:
        data = file.readlines()
    if not data[-1].endswith("\n"):
        data[-1] += "\n"   # add to avoid a messy append
    if include_N or added_sp2 == "amnt":
        data.append('no3'+'\n')

    # data.insert(1, 'fe2\t\n')
    with open(dst, 'w') as file:
        file.writelines(data)
    # remove duplicate minerals
    shf.remove_duplicates(dst)
        
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
    if not data[-1].endswith("\n"):
        data[-1] += "\n"   # add to avoid a messy append
    # add dust sp
    data.insert(1, added_sp+'\n')
    if added_sp2 in ['gbas','cc','cao','dlm','amnt', 'wls']: # then add this as well
        data.insert(1, added_sp2+'\n')
    if include_N and added_sp2 != "amnt":
        data.append('amnt'+'\n')
    if include_Al: 
        data.append(alphase+'\n')

    with open(dst, 'w') as file:
        file.writelines(data)
    # remove duplicate minerals
    shf.remove_duplicates(dst)

        
    filename = '/solutes.in'
    src = outdir + spinup_lab + filename
    dst = outdir + runname_lab  + filename
    with open(src, 'r') as file:
        data = file.readlines()
    if not data[-1].endswith("\n"):
        data[-1] += "\n"   # add to avoid a messy append
    if include_N or added_sp2 == "amnt":
        data.append('no3'+'\t\n')

    # data.insert(1, 'k\t\nna\t\nmg\t\n')
    with open(dst, 'w') as file:
        file.writelines(data)
    # remove duplicate minerals
    shf.remove_duplicates(dst)
        
    filename = '/kinspc.in'
    src = outdir + spinup_lab + filename
    dst = outdir + runname_lab  + filename
    with open(src, 'r') as file:
        data = file.readlines()
    if not data[-1].endswith("\n"):
        data[-1] += "\n"   # add to avoid a messy append
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
    
    
        
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ## --- setup for field run --- ##
    
    filename = '/frame.in'
    src = outdir + spinup_field + filename
    dst = outdir + runname_field  + filename

    with open(src, 'r') as file:
        data = file.readlines()
    data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(tau)
    if kwargs.get('mat') is not None:
        data[4] = '{:.8f}\ttemperature [oC]\n'.format(mat) 
    data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust) 
    data[6]     = '{:.8f}\tamounts of 2nd dusts [g/m2/yr]\n'.format(fdust2)
    data[7]     = '{:.8f}\tduration of dust application [yr]\n'.format(taudust)
    if kwargs.get('qrun') is not None:
        data[15] = '{:.8f}\tnet water flux [m/yr]\n'.format(mat)
    data[16]    = '{:.8f}\tradius of particles [m]\n'.format(dustrad)   # [tykukla added]
    data[18]    = '{}\n'.format(spinup_field)
    data[20]    = '{}\n'.format(runname_field)
    with open(dst, 'w') as file:
        file.writelines(data)
        
        
    ## --- run field run --- ##
    
    print(outdir+runname_field+'/scepter')
    os.system("chmod +x " + os.path.join(outdir,runname_field,'scepter'))  # grant permissions
    os.system(os.path.join(outdir,runname_field,'scepter'))


    ## --- getting data from field run --- ##
    
    # sps = ['g2','inrt',added_sp]
    sps = ['g2', added_sp]
    if include_Al:  sps.append( alphase )
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    
    aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
    # print(aqsps,btmconcs,dep)
    
    dic,dep = get_int_prof.get_ave_DIC_save(outdir,runname_field,dep_sample) # returning DIC in mol/ solid m3
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ## --- setup for lab run --- ##
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
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
        if sp not in oxide_ctnm_list:
            continue
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
    os.system("chmod +x " + os.path.join(outdir,runname_lab,'scepter'))  # grant permissions
    os.system(os.path.join(outdir,runname_lab,'scepter'))

    ## --- getting data from lab run --- ##
    
    # phint, ph_trend_lab = get_int_prof.get_ph_int_site_trend(outdir,runname_lab,dep_sample)
    # print(phint_field,phint)
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    # ... [TK] add a run=completed file for easier programmatic querying
    for runname in [runname_field,runname_lab]:
        dst = outdir + runname + where + 'completed.res'
        # Open the file in write mode to create it
        with open(dst, 'w') as f:
            pass  # pass does nothing, creating an empty file
        
    if use_local_storage:
        for runname in [runname_field,runname_lab]:
            src = outdir + runname 
            dst = outdir_src + runname 
            
            if not os.path.exists(dst): 
                shutil.copytree(src, dst)
            else:
                shutil.rmtree(dst)
                shutil.copytree(src, dst)

    
    # ... run postprocessing checks
    shf.run_complete_check(runname_field, 
                          runname_lab, 
                          outdir, 
                          target_duration=tau, 
                          include_duration_check=True, 
                          omit_saveSuff=True, 
                          omit_ipynb=True,
                         )

    # ... update the dust flux file so it's usable without needing input from frame.in
    shf.dustflx_calc(outdir, runname_field, fdust, fdust2, dustsp, dustsp_2nd)
    
    # UPDATE THE RUNNAME FOR NEXT ITERATION
    counter += 1
    runname_field_old = runname_field
    runname_lab_old = runname_lab
    # add to lists that are used for aws moving
    runname_field_list.append(runname_field)
    runname_lab_list.append(runname_lab)

    

# %%
# ---------------------------------------
# --- BUILD THE COMPOSITE DIRECTORY
print(expid) # (troubleshoot)
compDir_field, compDir_lab = cfxns.build_composite(expid, outdir)
runname_field_list.append(compDir_field)
runname_lab_list.append(compDir_lab)
# ---------------------------------------
# ... compute cdr-relevant fluxes
cflx.cflx_calc(outdir, compDir_field, [dustsp, dustsp_2nd])

# ... compute profile data
cflx.prof_postproc_save(outdir, runname_field, runname_lab, postproc_prof_list)

# ... move to aws if this option is turned on
# [nothing happens if aws_save != 'move' or 'copy']
outdir_postproc = shf.to_aws(aws_save, 
                            aws_bucket, 
                            outdir, 
                            runname_lab_list, 
                            runname_field_list)

# %%
