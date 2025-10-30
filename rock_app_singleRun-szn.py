# %% 
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
import subprocess

# --- read in helper functions from ew-workflows
from ew_workflows import scepter_helperFxns as shf
from ew_workflows import cflx_proc as cflx
# ---
# %% 

# -------------------------------------------------------------
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
tau = duration
added_sp = dustsp
added_sp2 = dustsp_2nd
spinid = spinrun
spinup_field = spinid+'_field'
spinup_lab   = spinid+'_lab'
expid = newrun_id
runname_field = expid+'_'+added_sp+'_field_tau'+str(tau).replace('.','p')
runname_lab   = expid+'_'+added_sp+'_lab_tau'+str(tau).replace('.','p')

if spinup_on: # then make sure the runnames end with "_field" or "_lab"
    print("spinup is on")
    runname_field = expid+'_'+added_sp+'_tau'+str(tau).replace('.','p')+'_field'
    runname_lab   = expid+'_'+added_sp+'_tau'+str(tau).replace('.','p')+'_lab'

fdust = dustrate
fdust2 = dustrate_2nd
dustrad = float(dustrad)/1e6  #  /1e6 converts micron to meters

datadir = os.path.join(modeldir, 'data/')

# --- files to delete if they are copied over
# (diagnostic files from the spinup run)
files_to_delete = ["check_logs.res", "check_results.res", "completed.res"]

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
    # save file denoting the variables used
    fn_dict_save = os.path.join(dst, "vars.res")
    shf.save_dict_to_text_file(combined_dict, fn_dict_save, delimiter='\t')
    # delete the diagnistic files from the spinup
    for file in files_to_delete:
        if os.path.exists(os.path.join(dst, file)):
            os.remove(os.path.join(dst, file))
    # make sure we have the right scepter run script
    shf.check_scepter_exec(scepter_exec_name, dst, modeldir)

# duplicate the climate files
if singlerun_seasonality==True:
    for runname in [runname_field,runname_lab]:
        # set source and destination
        src_clim = os.path.join(climatedir, climatefiles)
        dst_clim = os.path.join(outdir, runname)
        # copy over
        shf.copy_files(src_clim, dst_clim)

exename = scepter_exec_name
to = ' '
where = '/'

# ============ common input file modification wrt spinup: for field run
# --- UPDATE CEC IF ASKED ---------
if cec_update_from_spinup:
    for runname in [runname_field, runname_lab]:
        with open(os.path.join(outdir, runname_field, "cec.in"), "r") as f:
            lines = f.readlines()

        # Process each line
        modified_lines = []
        for line in lines:
            if re.match(r"^\S+\s+[-+]?\d*\.\d+", line):  # Check if line starts with a name followed by a number
                parts = line.split()
                parts[1] = str(cec)  # Modify first number
                parts[-1] = str(alpha)  # Modify last number
                modified_lines.append("\t".join(parts) + "\n")
            else:
                modified_lines.append(line)  # Keep header unchanged

        # print(modified_lines)

        # write back to the same file
        with open(os.path.join(outdir, runname, "cec.in"), "w") as f:
            f.writelines(modified_lines)
# -----------------------------------

filename = '/switches.in'
# ************************************************* 
# --- TROUBLESHOOT 
print("**************** TROUBLESHOOT *****************")
print(f"{spindir} ||| {spinup_field} ||| {filename}")
src = os.path.join(spindir, spinup_field, filename)
print(src)
print("***********************************************")
# ************************************************* 
dst = os.path.join(outdir, runname_field, filename)
with open(src, 'r') as file:
    data = file.readlines()
data[2] = '{:d}\tbio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken\n'.format(int(imix))
data[7] = 'true\trestart from a previous run\n'
if include_psd_bulk == True:
    data[-3] = 'true\tenabling PSD tracking\n'
elif include_psd_bulk == False:
    data[-3] = 'false\tenabling PSD tracking\n'
if include_psd_full == True:
    data[-2] = 'true\tenabling PSD tracking for individual solid species\n'
elif include_psd_full == False:
    data[-2] = 'false\tenabling PSD tracking for individual solid species\n'

if poro_iter_field == True:
    data[3] = 'true\tporosity  iteration\n'
elif poro_iter_field == False:
    data[3] = 'false\tporosity  iteration\n'
if poro_evol == True:
    data[-6] = 'true\tenabling porosity evolution\n'
elif poro_evol == False:
    data[-6] = 'false\tenabling porosity evolution\n'

if include_roughness_sa == True:
    data[8] = 'true\tinclude roughness in mineral surface area\n'
elif include_roughness_sa == False:
    data[8] = 'false\tinclude roughness in mineral surface area\n'
if sa_rule2 == True:
    data[-4] = 'true\tenabling SA evolution 2 (SA increases with porosity)\n'
elif sa_rule2 == False:
    data[-4] = 'false\tenabling SA evolution 2 (SA increases with porosity)\n'
if sa_rule1 == True:
    data[-5] = 'true\tenabling SA evolution 1 (SA decreases as porosity increases)\n'
elif sa_rule1 == False:
    data[-5] = 'false\tenabling SA evolution 1 (SA decreases as porosity increases)\n'

if singlerun_seasonality==True:
    data[-1] = 'true\tenabling seasonality\n'   # added per yoshi suggestion
elif singlerun_seasonality == False:
    data[-1] = 'false\tenabling full seasonality\n'
if cec_adsorption_on == True:
    data[11] = "true\tenabling adsorption for cation exchange\n"
elif cec_adsorption_on == False:
    data[11] = "false\tenabling adsorption for cation exchange\n"
    
with open(dst, 'w') as file:
    file.writelines(data)

# %% 

# tx = "test"
# if tx:
#     print("yay")
# elif not tx:
#     print("nay")

# %% 

# --- TROUBLESHOOT ------------ # 
print(src + "--------" + dst)
print(data)
time.sleep(5)


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
if added_sp2 == 'wls': dustsrc2 = os.path.join(modeldir, 'data', 'dust_wls.in')
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


# ------------------------------------------------------------------------------ # 
# *** added to test seasonality [TK - 4/10/24] *** #
# make_inputs.get_input_climate_temp(outdir = outdir,
#                                   runname = runname_field,
#                                   T_ave       = 15,
#                                   T_amp       = 0.3,
#                                   moist_ave   = 0.5,
#                                   moist_amp   = 0.3,
#                                   q_ave       = 0.5,
#                                   q_amp       = 0.3,
#                                   tau         = 1,
#                                    timeline =   np.linspace(0,1,12,endpoint=False)
#                                   )
# ------------------------------------------------------------------------------ # 

filename = '/slds.in'
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
filename = '/solutes.in'
src = os.path.join(spindir, spinup, filename)
dst = os.path.join(outdir, runname_field, filename)
with open(src, 'r') as file:
    data = file.readlines()
if not data[-1].endswith("\n"):
    data[-1] += "\n"   # add to avoid a messy append
if include_N or added_sp2 == "amnt":
    data.append('no3'+'\t\n')

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

filename = '/solutes.in'
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

filename = '/kinspc.in'
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
srcfile = os.path.join(modeldir, 'data', '2ndslds_def.in')
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

res_list = []
cnt = 0

# get ph of spin-up
phint = get_int_prof.get_ph_int_site(spindir,spinup_lab,dep_sample)
phint_field = get_int_prof.get_ph_int_site(spindir,spinup_field,dep_sample)

cnt += 1

    
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
## --- setup for field run --- ##

filename = '/frame.in'
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
    print(f"MAT: -------------------- {kwargs.get('mat')}")
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
with open(dst, 'w') as file:
    file.writelines(data)
    
    
## --- run field run --- ##

print(outdir+runname_field+'/' + exename)
os.system("chmod +x " + os.path.join(outdir,runname_field,exename))  # grant permissions
os.system(os.path.join(outdir,runname_field,exename))

## --- getting data from field run --- ##

phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)

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
    # print("THIS SPECIES ------" + sp) # TROUBLESHOOT
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
src = os.path.join(spindir, spinup_lab, filename)
dst = os.path.join(outdir, runname_lab, filename)

with open(src, 'r') as file:
    data = file.readlines()
data[1]     = '{:.8f}\ttotal depth of weathering profile [m]\n'.format(ztot_lab)
data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(ttot_lab)
data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust_lab)
data[16]    = '{:.8f}\tradius of particles [m]\n'.format(dustrad)   # [tykukla added]
data[10]     = '{:.8f}\tinitial porosity\n'.format(poro_lab)
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

print(outdir+runname_lab+'/' +exename)
os.system("chmod +x " + os.path.join(outdir,runname_lab,exename))  # grant permissions
os.system(os.path.join(outdir,runname_lab,exename))

## --- getting data from lab run --- ##

phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
print(phint_field,phint)

## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

print(' ')
print(' ')
print(' ')
print(' ')
# print("******* {:d}'s iteration, error = {:.8f}".format(cnt,error))
print(' ')
print(' ')
print(' ')
print(' ')

res_list.append([cnt, phint_field, phint, fdust ])

cnt += 1


time.sleep(5)

for runname in [runname_field,runname_lab]:
    np.savetxt(outdir + runname + where + 'iteration_tmp.res',np.array(res_list))



name_list = [
    'iter.'
    ,'porewater_pH[-]'
    ,'soil_pHw[-]'
    ,'basalt[g/m2/yr]'
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

# ... compute cdr-relevant fluxes
multi_sp_dict = {dustsp: multi_sp_feedstock, dustsp_2nd: multi_sp_feedstock_2nd}
cflx.cflx_calc(outdir, runname_field, [dustsp, dustsp_2nd], multi_sp_dict)

# ... compute profile data
cflx.prof_postproc_save(outdir, runname_field, runname_lab, postproc_prof_list)

# ... move to aws if this option is turned on
# [nothing happens if aws_save != 'move' or 'copy']
shf.to_aws(aws_save, 
           aws_bucket, 
           outdir, 
           runname_lab, 
           runname_field)



# %%
