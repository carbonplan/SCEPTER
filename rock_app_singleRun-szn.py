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

# --- read in helper functions from aglime-swap-cdr
# add aglime-swap-cdr dir to path # [UPDATE FOR YOUR MACHINE]
sys.path.append(os.path.abspath('/home/tykukla/aglime-swap-cdr/scepter/setup'))
# import module
import scepter_helperFxns as shf
# ---


# -------------------------------------------------------------
# -------------------------------------------------------------
# --- set default and system args
sys_args = shf.parse_arguments(sys.argv)   # parse system args
import_dict = sys_args['default_dict'] # set the dictionary to use from system args
def_args = getattr(defaults.dict_singlerun, import_dict)  # get dict attribute

# set global variables
combined_dict = shf.set_vars(def_args, sys_args)  # default unless defined in sys_args
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
fdust = dustrate
fdust2 = dustrate_2nd
dustrad = float(dustrad)/1e6  #  /1e6 converts micron to meters

datadir = os.path.join(modeldir, 'data/')

# duplicate directories from spinups
for (runname,spinup) in [(runname_field,spinup_field),(runname_lab,spinup_lab)]:

    src = outdir + spinup
    dst = outdir + runname

    if not os.path.exists(dst): 
        shutil.copytree(src, dst)
    else:
        shutil.rmtree(dst)
        shutil.copytree(src, dst)
    # save file denoting the variables used
    fn_dict_save = os.path.join(dst, "vars.res")
    shf.save_dict_to_text_file(combined_dict, fn_dict_save, delimiter='\t')

# duplicate the climate files
if singlerun_seasonality==True:
    for runname in [runname_field,runname_lab]:
        # set source and destination
        src_clim = os.path.join(climatedir, climatefiles)
        dst_clim = os.path.join(outdir, runname)
        # copy over
        shf.copy_files(src_clim, dst_clim)

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


# --- TROUBLESHOOT ------------ # 
print(src + "--------" + dst)
print(data)
time.sleep(5)


# --- PRIMARY DUST FILE
if added_sp == 'gbas': dustsrc = os.path.join(modeldir, 'data', 'dust_gbasalt.in')
if added_sp == 'cc': dustsrc = os.path.join(modeldir, 'data', 'dust_lime.in')
if added_sp == 'cao': dustsrc = os.path.join(modeldir, 'data', 'dust_cao.in')
if added_sp == 'dlm': dustsrc = os.path.join(modeldir, 'data', 'dust_dlm.in')
dustdst = 'dust.in'

os.system('cp ' + dustsrc + to + outdir + runname_field + where + dustdst) 


# --- SECONDARY DUST FILE
if added_sp2 == "amnt": dustsrc2 = os.path.join(modeldir, 'data', 'dust_fert.in')
if added_sp2 == 'gbas': dustsrc2 = os.path.join(modeldir, 'data', 'dust_gbasalt.in')
if added_sp2 == 'cc': dustsrc2 = os.path.join(modeldir, 'data', 'dust_lime.in')
if added_sp2 == 'cao': dustsrc2 = os.path.join(modeldir, 'data', 'dust_cao.in')
if added_sp2 == 'dlm': dustsrc2 = os.path.join(modeldir, 'data', 'dust_dlm.in')  
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
src = outdir + spinup_field + filename
dst = outdir + runname_field  + filename
with open(src, 'r') as file:
    data = file.readlines()

data.insert(1, added_sp+'\t\n')
if include_N:
    data.append('amnt'+'\t\n')
if include_Al: 
    data.append(alphase+'\t\n')
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
if include_N:
    data.append('no3'+'\t\n')

# data.insert(1, 'fe2\t\n')
with open(dst, 'w') as file:
    file.writelines(data)
    
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
if include_N:
    data.append('amnt'+'\t\n')
if include_Al: 
    data.append(alphase+'\t\n')
with open(dst, 'w') as file:
    file.writelines(data)
# remove duplicate minerals
shf.remove_duplicates(dst)

filename = '/solutes.in'
src = outdir + spinup_lab + filename
dst = outdir + runname_lab  + filename
with open(src, 'r') as file:
    data = file.readlines()
if include_N:
    data.append('no3'+'\t\n')
# data.insert(1, 'k\t\nna\t\nmg\t\n')
with open(dst, 'w') as file:
    file.writelines(data)

filename = '/kinspc.in'
src = outdir + spinup_lab + filename
dst = outdir + runname_lab  + filename
with open(src, 'r') as file:
    data = file.readlines()
    
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
phint = get_int_prof.get_ph_int_site(outdir,spinup_lab,dep_sample)
phint_field = get_int_prof.get_ph_int_site(outdir,spinup_field,dep_sample)

cnt += 1

    
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
## --- setup for field run --- ##

filename = '/frame.in'
src = outdir + spinup_field + filename
dst = outdir + runname_field  + filename

with open(src, 'r') as file:
    data = file.readlines()
data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(tau)
data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust) 
data[6]     = '{:.8f}\tamounts of 2nd dusts [g/m2/yr]\n'.format(fdust2)
data[7]     = '{:.8f}\tduration of dust application [yr]\n'.format(taudust)
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

phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
sps = ['g2',added_sp] # ['g2','inrt',added_sp]    # TK -- come back to these so all solid species are considered? is it necessary?
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
src = outdir + spinup_lab + filename
dst = outdir + runname_lab  + filename

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

print(outdir+runname_lab+'/scepter')
os.system("chmod +x " + os.path.join(outdir,runname_lab,'scepter'))  # grant permissions
os.system(os.path.join(outdir,runname_lab,'scepter'))

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


# ... move to aws if this option is turned on
# [nothing happens if aws_save != 'move' or 'copy']
shf.to_aws(aws_save, 
           aws_bucket, 
           outdir, 
           runname_lab, 
           runname_field)


