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
import fsspec
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
from ew_workflows import scepter_helperFxns as shf
from ew_workflows import build_composite_multiyear as cfxns
from ew_workflows import cflx_proc as cflx
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
        spindir = outdir # set the spinup dir to outdir for later iterations
        # fdust = next_dustrate  # set fdust to some ideally smaller value for the successive iterations
    
    # set runnames
    runname_field   = f"{expid}_startyear-{str(tstep).replace('.','p')}_iter-{int(counter)}_field"
    runname_lab     = f"{expid}_startyear-{str(tstep).replace('.','p')}_iter-{int(counter)}_lab"

    # --- move spinup files to local if on s3
    if spindir.startswith("s3://"):
        spindir = shf.copy_spinup_if_on_s3( # update spindir to local
                spindir=spindir, 
                spinup_field=spinup_field, 
                spinup_lab=spinup_lab,
                localdir=outdir,
        )
    # ---

    # duplicate directories from spinups
    for (runname,spinup) in [(runname_field,spinup_field),(runname_lab,spinup_lab)]:
        
        src = os.path.join(spindir, spinup)
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
        # make sure we have the right scepter run script
        shf.check_scepter_exec(scepter_exec_name, dst, modeldir)
    
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
    exename = scepter_exec_name
    to = ' '
    where = '/'
    
    # ============ common input file modification wrt spinup: for field run
    filename = 'switches.in'
    src = os.path.join(spindir, spinup_field, filename)
    dst = os.path.join(outdir, runname_field, filename)
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
    
    # === modifications required for v1.0.2 ----------------------
    if exename in v102_exelist:
        data[5] = '{:d}\tdisplay results at runtime: 0-- none, 1-- only reporting time, 2-- every time iteration, if not defined 1 is taken\n'.format(int(1))
        data[6] = '{:d}\treport files: 0-- basics, 1-- +saturation time series\n'.format(int(1))


    with open(dst, 'w') as file:
        file.writelines(data)
    

    # --- write Dust_temp.in (for v1.0.2 seasonal runs)
    if exename in v102_exelist and singlerun_seasonality:
        shf.create_dust_input(
            outdir = outdir,
            runname = runname_field,
            dustname1 = added_sp,
            dustname2 = added_sp2,
            t_add = duststart,
            output_filename  = "Dust_temp.in",
            dryrun = False
        )
    # --------------------------------------------------


    # --- PRIMARY DUST FILE
    multi_sp_feedstock = False
    if added_sp == "amnt": dustsrc = os.path.join(modeldir, 'data', 'dust_fert.in')
    if added_sp == 'gbas': dustsrc = os.path.join(modeldir, 'data', 'dust_gbasalt.in')
    if added_sp == 'cc': dustsrc = os.path.join(modeldir, 'data', 'dust_lime.in')
    if added_sp == 'cao': dustsrc = os.path.join(modeldir, 'data', 'dust_cao.in')
    if added_sp == 'dlm': dustsrc = os.path.join(modeldir, 'data', 'dust_dlm.in')
    if added_sp == 'wls': dustsrc = os.path.join(modeldir, 'data', 'dust_wls.in')
    if added_sp == 'fo': dustsrc = os.path.join(modeldir, 'data', 'dust_fo.in')
    if added_sp == 'baek23': 
        dustsrc = os.path.join(modeldir, 'data', 'dust_def.in')
        multi_sp_feedstock = True
    if added_sp == 'baek23nodp': 
        dustsrc = os.path.join(modeldir, 'data', 'dust_def_noDP.in')
        multi_sp_feedstock = True
    if added_sp == 'bridge': # blue ridge basalt
        dustsrc = os.path.join(modeldir, 'data', 'dust_BlueRidge.in')
        multi_sp_feedstock = True
    dustdst = 'dust.in'
    
    os.system('cp ' + dustsrc + to + outdir + runname_field + where + dustdst) 

    # --- SECONDARY DUST FILE
    multi_sp_feedstock_2nd = False
    if added_sp2 == "amnt": dustsrc2 = os.path.join(modeldir, 'data', 'dust_fert.in')
    if added_sp2 == 'gbas': dustsrc2 = os.path.join(modeldir, 'data', 'dust_gbasalt.in')
    if added_sp2 == 'cc': dustsrc2 = os.path.join(modeldir, 'data', 'dust_lime.in')
    if added_sp2 == 'cao': dustsrc2 = os.path.join(modeldir, 'data', 'dust_cao.in')
    if added_sp2 == 'dlm': dustsrc2 = os.path.join(modeldir, 'data', 'dust_dlm.in')  
    if added_sp2 == 'wls': dustsrc = os.path.join(modeldir, 'data', 'dust_wls.in')
    if added_sp2 == 'fo': dustsrc2 = os.path.join(modeldir, 'data', 'dust_fo.in')
    if added_sp2 == 'baek23': 
        dustsrc2 = os.path.join(modeldir, 'data', 'dust_def.in')
        multi_sp_feedstock_2nd = True
    if added_sp2 == 'baek23nodp': 
        dustsrc2 = os.path.join(modeldir, 'data', 'dust_def_noDP.in')
        multi_sp_feedstock_2nd = True
    if added_sp2 == 'bridge': # blue ridge basalt
        dustsrc2 = os.path.join(modeldir, 'data', 'dust_BlueRidge.in')
        multi_sp_feedstock_2nd = True
    dustdst2 = 'dust_2nd.in'
    
    os.system('cp ' + dustsrc2 + to + outdir + runname_field + where + dustdst2) 
    # --- 
    
    filename = 'slds.in'
    src = os.path.join(spindir, spinup_field, filename)
    dst = os.path.join(outdir, runname_field, filename)
    with open(src, 'r') as file:
        data = file.readlines()
    if not data[-1].endswith("\n"):
        data[-1] += "\n"   # add to avoid a messy append
    # add dust species if name is mineral
    if not multi_sp_feedstock:
        data.insert(1, added_sp+'\n')
    else: # otherwise get dustsp from the dust.in file
        data = shf.add_dustsp_to_sld(data, dustdst, outdir, runname_field)
        
    # add second species if needed
    if not multi_sp_feedstock_2nd: # then add this as well
        data.insert(1, added_sp2+'\n')
    else: # otherwise get dustsp from the dust.in file
        data = shf.add_dustsp_to_sld(data, dustdst2, outdir, runname_field)
    if include_N and added_sp2 != "amnt":
        data.append('amnt'+'\n')
    if include_Al: 
        data.append(alphase+'\n')

    with open(dst, 'w') as file:
        file.writelines(data)
    # remove duplicate minerals
    shf.remove_duplicates(dst)

        
    # ============ adding Fe(II) as tracer and its oxidation =================
    filename = 'solutes.in'
    src = os.path.join(spindir, spinup, filename)
    dst = os.path.join(outdir, runname_field, filename)
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
        
    # filename = 'gases.in'
    # src = outdir + spinup + filename
    # dst = outdir + runname  + filename
    # with open(src, 'r') as file:
        # data = file.readlines()
    
    # data.insert(1, 'po2\t\n')
    # with open(dst, 'w') as file:
        # file.writelines(data)
        
    # filename = 'extrxns.in'
    # src = outdir + spinup + filename
    # dst = outdir + runname  + filename
    # with open(src, 'r') as file:
        # data = file.readlines()
    
    # data.insert(1, 'fe2o2\t\n')
    # with open(dst, 'w') as file:
        # file.writelines(data)
        
    # ============ common input file modification wrt spinup: for lab run ============
    
    filename = 'slds.in'
    src = os.path.join(spindir, spinup_lab, filename)
    dst = os.path.join(outdir, runname_lab, filename)
    with open(src, 'r') as file:
        data = file.readlines()
    if not data[-1].endswith("\n"):
        data[-1] += "\n"   # add to avoid a messy append
    # add dust species if name is mineral
    if not multi_sp_feedstock:
        data.insert(1, added_sp+'\n')
    else: # otherwise get dustsp from the dust.in file
        data = shf.add_dustsp_to_sld(data, dustdst, outdir, runname_lab)
        
    # add second species if needed
    if not multi_sp_feedstock_2nd: # then add this as well
            data.insert(1, added_sp2+'\n')
    else: # otherwise get dustsp from the dust.in file
        data = shf.add_dustsp_to_sld(data, dustdst2, outdir, runname_lab)
    if include_N and added_sp2 != "amnt":
        data.append('amnt'+'\n')
    if include_Al: 
        data.append(alphase+'\n')

    with open(dst, 'w') as file:
        file.writelines(data)
    # remove duplicate minerals
    shf.remove_duplicates(dst)

        
    filename = 'solutes.in'
    src = os.path.join(spindir, spinup_lab, filename)
    dst = os.path.join(outdir, runname_lab, filename)
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
        
    filename = 'kinspc.in'
    src = os.path.join(spindir, spinup_lab, filename)
    dst = os.path.join(outdir, runname_lab, filename)
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

    # --- particle size distribution setup
    if include_psd_bulk or include_psd_full:
        # define the source file 
        if use_psdrain_datfile:
            srcfile = os.path.join(modeldir, 'data', psdrain_datfile)
        else: 
            srcfile = None  # this means we make one from scratch using sld_varlist
        
        filename = 'psdrain.in'
        sld_varlist = [[psdrain_meanRad, psdrain_log10_sd, psdrain_wt]]
        # make (or copy) psdrain.in
        make_inputs.get_input_sld_properties(
            outdir = outdir
            ,runname = runname_field
            ,filename = filename 
            ,srcfile = srcfile
            ,sld_varlist = sld_varlist,
        )
    # -------------------------------------
    # compile 
    # os.system('make')
    # os.system('make --file=makefile_test')
    # for runname in [runname_field,runname_lab]:
        # if not os.path.exists( outdir + runname) : 
            # os.system('mkdir -p ' + outdir + runname)
            # os.system('cp ' + exename_src + to + outdir + runname + where + exename)
    
    
        
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ## --- setup for field run --- ##
    
    filename = 'frame.in'
    src = os.path.join(spindir, spinup_field, filename)
    dst = os.path.join(outdir, runname_field, filename)

    with open(src, 'r') as file:
        data = file.readlines()
    if kwargs.get('ztot_field') is not None and not np.isnan(kwargs.get('ztot_field')):
        print(f"z_field: -------------------- {kwargs.get('ztot_field')}")
        data[1] = '{:.8f}\ttotal depth of weathering profile [m]'.format(ztot_field)
    if kwargs.get('nz') is not None and not np.isnan(kwargs.get('nz')):
        print(f"nz: -------------------- {kwargs.get('nz')}")
        data[2] = '{:.8f}\ttotal depth of weathering profile [m]'.format(nz)
    data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(tau)
    if kwargs.get('mat') is not None and not np.isnan(kwargs.get('mat')):
        data[4] = '{:.8f}\ttemperature [oC]\n'.format(mat) 
    data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust) 
    data[6]     = '{:.8f}\tamounts of 2nd dusts [g/m2/yr]\n'.format(fdust2)
    data[7]     = '{:.8f}\tduration of dust application [yr]\n'.format(taudust)
    if kwargs.get('soilmoisture_surf') is not None:
        data[11] = '{:.8f}\twater saturation at the surface of profile\n'.format(soilmoisture_surf)
    if kwargs.get('dust_mixdep') is not None:
        data[13] = '{:.8f}\tdepth of mixed layer for dust [m]\n'.format(dust_mixdep)
    if kwargs.get('qrun') is not None:
        data[15] = '{:.8f}\tnet water flux [m/yr]\n'.format(qrun)
    data[16]    = '{:.8f}\tradius of particles [m]\n'.format(dustrad)   # [tykukla added]
    data[18]    = '{}\n'.format(spinup_field)
    data[20]    = '{}\n'.format(runname_field)

    # --- if v102 then use full path instead of just spinname
    if exename in v102_exelist:
        data[18]    = '{}\n'.format(os.path.join(outdir, spinup_field))
    # ---

    with open(dst, 'w') as file:
        file.writelines(data)
        
        
    ## --- run field run --- ##
    
    print(outdir+runname_field+'/' + exename)
    os.system("chmod +x " + os.path.join(outdir,runname_field,exename))  # grant permissions
    os.system(os.path.join(outdir,runname_field,exename))


    ## --- getting data from field run --- ##
    
    if not multi_sp_feedstock:
        sps = ['g2',added_sp] # ['g2','inrt',added_sp]    # TK -- come back to these so all solid species are considered? is it necessary?
    else: # otherwise get dustsp from the dust.in file
        slds_path = os.path.join(outdir, runname_field, 'slds.in')
        with open(slds_path, "r") as f:
            sld_lines = f.readlines()
        # undo formatting
        sps = [line.strip() for line in sld_lines[1:]]
        
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
    
    filename = 'frame.in'
    src = os.path.join(spindir, spinup_lab, filename)
    dst = os.path.join(outdir, runname_lab, filename)

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
    

    filename = 'switches.in'
    src = os.path.join(spindir, spinup_lab, filename)
    dst = os.path.join(outdir, runname_lab, filename)

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
    
    # === modifications required for v1.0.2 ----------------------
    if exename in v102_exelist:
        data[5] = '{:d}\tdisplay results at runtime: 0-- none, 1-- only reporting time, 2-- every time iteration, if not defined 1 is taken\n'.format(int(1))
        data[6] = '{:d}\treport files: 0-- basics, 1-- +saturation time series\n'.format(int(1))


    with open(dst, 'w') as file:
        file.writelines(data)
    # break
    
    # --- write Dust_temp.in (for v1.0.2 seasonal runs)
    if exename in v102_exelist and singlerun_seasonality:
        shf.create_dust_input(
            outdir = outdir,
            runname = runname_field,
            dustname1 = added_sp,
            dustname2 = added_sp2,
            t_add = duststart,
            output_filename  = "Dust_temp.in",
            dryrun = False
        )
    # --------------------------------------------------


    ## --- run lab run --- ##
    
    print(outdir+runname_lab+'/'+exename)
    os.system("chmod +x " + os.path.join(outdir,runname_lab,exename))  # grant permissions
    os.system(os.path.join(outdir,runname_lab,exename))

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

# ... run postprocessing checks
shf.run_complete_check(compDir_field, 
                      compDir_lab, 
                      outdir, 
                      target_duration=tau, 
                      include_duration_check=False, 
                      omit_saveSuff=True, 
                      omit_ipynb=True,
                     )

# %% 
# ... compute cdr-relevant fluxes
multi_sp_dict = {dustsp: multi_sp_feedstock, dustsp_2nd: multi_sp_feedstock_2nd}
cflx.cflx_calc(outdir, compDir_field, [dustsp, dustsp_2nd], multi_sp_dict)
# %% 
# ... compute profile data
cflx.prof_postproc_save(outdir, compDir_field, compDir_lab, postproc_prof_list, multi_iter=True)
# %% 
# ... move to aws if this option is turned on
# [nothing happens if aws_save != 'move' or 'copy']
outdir_postproc = shf.to_aws(aws_save, 
                            aws_bucket, 
                            outdir, 
                            runname_lab_list, 
                            runname_field_list)

# %%
