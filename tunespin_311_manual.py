import os
import numpy as np
import get_int_prof
import make_inputs
import get_inputs
import json
import time
import sys
import re
from datetime import datetime
import shutil
import subprocess
import defaults.dict_singlerun

datestr = datetime.today().strftime('%Y-%m-%d')




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
    # save dict
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
# -------------------------------------------------------------
# -------------------------------------------------------------

print(sys.argv)

# --- set default and system args (TK)
sys_args = parse_arguments(sys.argv)   # parse system args
import_dict = "s311_dictionary" # sys_args['default_dict'] # set the dictionary to use from system args
def_args = getattr(defaults.dict_singlerun, import_dict)  # get dict attribute

# set global variables
# (TK: currently, combined_dict doesn't get saved in this version)
combined_dict = set_vars(def_args, sys_args)  # default unless defined in sys_args

water_frac = water_frac_tunespin
targetpH = tph
acidsat = tec
targetOM = tsom
poro_input = poro
moist_input = soilmoisture
sat_input = moist_input/poro
kh = 10**(3.5 + alpha/2.)

# --- TUNED INPUTS (TK added)
# omrain_field= 1 * 900
ca = 0 #  1 * 5e-3
kh = 1 * kh
# ---

runname_guess = initguess

datadir = os.path.join(modeldir, 'data/')

if liming: ca = 10

outdir = outdir
simid = spinname
runname_field   = simid+'_spintuneupManual_field'
runname_lab     = simid+'_spintuneupManual_lab'




# compile 
exename = 'scepter'
exename_src = 'scepter'
# exename_src = 'scepter_test'
to = ' '
where = '/'

# get current dur
mycwd = os.getcwd()
# os.system('make')
os.chdir(modeldir)  # change to model directory where makefile resides
os.system('make')




prev_iter_exist = False
if not os.path.exists( outdir + runname_field) : 
    os.system('mkdir -p ' + outdir + runname_field)
else:
    if os.path.exists( outdir + runname_field + where + 'iteration_tmp.res'):
        iter_prev = np.loadtxt(outdir + runname_field + where + 'iteration_tmp.res',skiprows=1)
        shutil.copy(outdir + runname_field + where + 'iteration_tmp.res'
            ,outdir + runname_field + where + 'iteration_tmp_SAVE_'+datestr+'.res')
        # prev_iter_exist = True
if not os.path.exists( outdir + runname_lab) : os.system('mkdir -p ' + outdir + runname_lab)
os.system('cp ' + exename_src + to + outdir + runname_field + where + exename)
os.system('cp ' + exename_src + to + outdir + runname_lab + where + exename)
# change back
os.chdir(mycwd)





# --- FRAME.IN SELECTIONS
# 
ztot=0.5
ztot_field=0.5
# ztot_field=0.3
ztot_lab=0.5
ztot_lab=0.05
nz=30
ttot_field=100000
ttot_lab=1000
ttot_lab=100
temp_field=mat
temp_lab=25
fdust_field=0
fdust_lab=0
fdust2=0
taudust_field=0
taudust_lab=0.01
omrain_lab=0
poro_field=poro_input
poro_lab=0.928391508
moistsrf_field=sat_input
moistsrf_lab=1.0
# zdust_field=0.18
zdust_field=0.25
zdust_lab=0.15
w_field=100e-5
w_field=erosion
w_lab=0
q_field=1200e-3
q_field=qrun
q_lab=0
p=10e-6
nstep=10
rstrt='self'
runid_field=runname_field
runid_lab=runname_lab

# N_rain = 0  # gN/m2/yr
# N_rain = 8.406375  # gN/m2/yr ( <---> 75 lbs/acre/year)
# N_rain = 24.6587   # gN/m2/yr ( <---> 220 lbs/acre/year)
N_rain = nitrif   # gN/m2/yr ( <---> 220 lbs/acre/year)
N_rain = N_rain/14  # mol N/m2/yr
N_rain = N_rain*80  # g NH4NO3/m2/yr
N_rain = N_rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
N_rain = N_rain/omrain_field  # normalize against OC rain 

CC_rain = ca # gCaCO3/m2/yr
CC_rain = CC_rain/omrain_field




# --- WRITE THE FRAME.IN FILE 
# 
make_inputs.get_input_frame(
    outdir=outdir
    ,runname=runname_field
    ,ztot=ztot_field
    ,nz=nz
    ,ttot=ttot_field
    ,temp=temp_field
    ,fdust=fdust_field
    ,fdust2=fdust2
    ,taudust=taudust_field
    ,omrain=omrain_field
    ,zom=zom
    ,poro=poro_field
    ,moistsrf=moistsrf_field
    ,zwater=zwater
    ,zdust=zdust_field
    ,w=w_field
    ,q=q_field
    ,p=p
    ,nstep=nstep
    ,rstrt=rstrt
    ,runid=runname_field
    )




# --- SWITCHES.IN SELECTIONS
# 
display='true'
disp_lim='true'
restart ='false'
act_ON ='false'
if activity_on : act_ON ='true'
dt_fix='false'
if cec_adsorption_on == True:
    cec_on="true"
else:
    cec_on="false"
dz_fix='true'
close_aq_field='false'
poro_evol='true'
sa_evol_1 ='true'
sa_evol_2='false'
season='false'
if include_roughness_sa == True:
    rough_field      ='true'
else:
    rough_field      ='false'

w_scheme_lab=0
mix_scheme_lab=0 
poro_iter_lab='true' 
rough_lab      ='false'
close_aq_lab='true'
psd_bulk_lab='false'
psd_full_lab ='false'




# --- WRITE THE SWITCHES.IN FILE 
# 
make_inputs.get_input_switches(
    outdir=outdir
    ,runname=runname_field
    ,w_scheme=w_scheme_field
    ,mix_scheme=mix_scheme_field
    ,poro_iter=poro_iter_field
    ,sldmin_lim=sldmin_lim 
    ,display=display
    ,disp_lim=disp_lim
    ,restart=restart 
    ,rough=rough_field
    ,act_ON=act_ON 
    ,dt_fix=dt_fix
    ,cec_on=cec_on
    ,dz_fix=dz_fix
    ,close_aq=close_aq_field
    ,poro_evol=poro_evol
    ,sa_evol_1=sa_evol_1 
    ,sa_evol_2=sa_evol_2
    ,psd_bulk=psd_bulk_field
    ,psd_full=psd_full_field
    ,season=season
    )




# --- WRITE INPUT SPECIATION:
# 
#     extrxns.in
#     gases.in
#     solutes.in
#     slds.in
#     
sld_list_field=['amsi', 'cc', 'dlm', 'tm', 'ka', 'dp', 'qtz', 'ill', 'cabd', 'kfs', 'ab', 'g2']
if liming:
    sld_list_field.append(limesp)
aq_list_field = ['mg', 'si', 'na', 'ca', 'al', 'fe3', 'k']
if include_N: 
    sld_list_field.append('amnt')
    aq_list_field.append('no3')
if include_Al: 
    sld_list_field.append(alphase)
    aq_list_field.append('al')
if add_secondary:                # [tk added to track secondary precipitates]
    sld_list_field += sld_track
    
gas_list_field = ['pco2', 'po2']
exrxn_list_field = []
# exrxn_list_field = ['g2ca']
make_inputs.get_input_tracers(
    outdir=outdir
    ,runname=runname_field
    ,sld_list = sld_list_field
    ,aq_list = aq_list_field
    ,gas_list = gas_list_field
    ,exrxn_list = exrxn_list_field
    )




# --- WRITE BC CONDITIONS
# 
#     parentrock.in
#     atm.in
#     rain.in
# 

# bring in parent rock file [tk added]
if (spinup_parentrock_file != None) & ("json" in spinup_parentrock_file):
    pr_file = os.path.join(modeldir, 'data', spinup_parentrock_file)
    with open(pr_file, 'r') as f:
        pr_list_field = json.load(f)
else:
    pr_list_field = [('inrt',1.0)]
    
rain_list_field = []  # [('ca',ca)]
if liming: rain_list_field = []
atm_list_field = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
if include_Al:  pr_list_field.append( (alphase,1e-2) )
make_inputs.get_input_tracer_bounds(
    outdir=outdir
    ,runname=runname_field
    ,pr_list = pr_list_field
    ,rain_list=rain_list_field
    ,atm_list=atm_list_field
    )




# --- WRITE DUST INPUTS
# 
#     dust.in
#     cec.in
#     psdrain.in
#     2ndslds.in
#     OM_rain.in
#     

filename = 'dust.in'
srcfile = os.path.join(datadir, 'dust_gbasalt.in')
sld_varlist =[]
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    # ,srcfile = srcfile
    ,sld_varlist=sld_varlist
    )
    
filename = 'cec.in'
sld_varlist = [('inrt',4,4.1,alpha) ,('g2',4,4.1,alpha) ] 
if include_Al:  sld_varlist.append( (alphase,4,4.1,alpha)  )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = 'OM_rain.in'
sld_varlist = [ ('g2',1) ] 
if include_N: sld_varlist.append( ('amnt',N_rain) ) 
if liming: sld_varlist.append( (limesp,CC_rain) ) 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
    )

# [TK ADDED THIS -- for some reason kinsp was only being added to in the lab run..]
filename = 'kinspc.in'
sld_varlist = [ ('g2',g2_turnoverTime) ]    # value used in v0.9 example 3-2
# sld_varlist = [ ('g2',1600) ] 
if include_Al:  sld_varlist.append(  (alphase,0)  )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
    )


filename = '2ndslds.in'
srcfile = os.path.join(datadir, '2ndslds_def.in')
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,srcfile = srcfile
    )
    
filename = 'psdrain.in'
srcfile = os.path.join(datadir, 'psdrain_320um.in')
sld_varlist = [ (5e-6,0.2,1), (5e-6,0.2,1), (5e-6,0.2,1), (5e-6,0.2,1), ] 
# sld_varlist = [ (1e-6,0.2,1), (2e-6,0.2,1), (5e-6,0.2,1), (20e-6,0.2,1), ] 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    # ,sld_varlist=sld_varlist
    ,srcfile = srcfile
    )
    
if liming:
    # to simplify/speed-up calculation
    filename = 'nopsd.in'
    sld_varlist = [ (limesp,) ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    # to simplify/speed-up calculation
    filename = 'sa.in'
    sld_varlist = [ (limesp,1e-5) ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )





# ============= buffer exp. setup ============= 
make_inputs.get_input_frame(
    outdir=outdir
    ,runname=runname_lab
    ,ztot=ztot_lab
    ,nz=nz
    ,ttot=ttot_lab
    ,temp=temp_lab
    ,fdust=fdust_lab
    ,fdust2=fdust2
    ,taudust=taudust_lab
    ,omrain=omrain_lab
    ,zom=zom
    ,poro=poro_lab
    ,moistsrf=moistsrf_lab
    ,zwater=zwater
    ,zdust=zdust_lab
    ,w=w_lab
    ,q=q_lab
    ,p=p
    ,nstep=nstep
    ,rstrt=rstrt
    ,runid=runname_lab
    )
    
make_inputs.get_input_switches(
    outdir=outdir
    ,runname=runname_lab
    ,w_scheme=w_scheme_lab
    ,mix_scheme=mix_scheme_lab
    ,poro_iter=poro_iter_lab
    ,sldmin_lim=sldmin_lim 
    ,display=display
    ,disp_lim=disp_lim
    ,restart=restart 
    ,rough=rough_lab
    ,act_ON=act_ON 
    ,dt_fix=dt_fix
    ,cec_on=cec_on
    ,dz_fix=dz_fix
    ,close_aq=close_aq_lab
    ,poro_evol=poro_evol
    ,sa_evol_1=sa_evol_1 
    ,sa_evol_2=sa_evol_2
    ,psd_bulk=psd_bulk_lab
    ,psd_full=psd_full_lab
    ,season=season
    )




sld_list_lab = ['inrt','g2','cao','mgo','na2o','k2o']
aq_list_lab = ['ca','mg','na','k'] 
gas_list_lab = [] # (added 3.23.2023)
if include_N: 
    sld_list_lab.append('amnt')
    aq_list_lab.append('no3')
if include_Al: 
    sld_list_lab.append(alphase)
    aq_list_lab.append('al')
if use_CaCl2:
    sld_list_lab.append('cacl2')
    aq_list_lab.append('cl')
if include_DIC: # (added 3.22.2023)
    sld_list_lab.append('g1')
    gas_list_lab.append('pco2')
make_inputs.get_input_tracers(
    outdir=outdir
    ,runname=runname_lab
    ,sld_list = sld_list_lab
    ,aq_list = aq_list_lab
    ,gas_list = gas_list_lab # (added 3.23.2023)
    # ,exrxn_list = exrxn_list
    )
    
# pr_list = [('inrt',1.0)]
# bring in parent rock file [tk added]
if (spinup_parentrock_file != None) & ("json" in spinup_parentrock_file):
    pr_file = os.path.join(modeldir, 'data', spinup_parentrock_file)
    with open(pr_file, 'r') as f:
        pr_list_field = json.load(f)
else:
    pr_list_field = [('inrt',1.0)]
pr_list_lab = [('inrt',0.95),('g2',0.05)]
atm_list_lab = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
if include_DIC: atm_list_lab = [('pco2',1e-20),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)] # (added 3.22.2023)
make_inputs.get_input_tracer_bounds(
    outdir=outdir
    ,runname=runname_lab
    ,pr_list = pr_list_lab
    # ,rain_list=rain_list
    ,atm_list=atm_list_lab
    )
    
filename = 'dust.in'
srcfile = os.path.join(datadir, 'dust_cao.in')
sld_varlist =[]
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,srcfile = srcfile
    # ,sld_varlist=sld_varlist
    )
    
filename = 'cec.in'
# sld_varlist = [('inrt',4,4.1)  ] 
sld_varlist = [('inrt',4,4.1,alpha) ,('g2',4,4.1,alpha) ] 
if include_Al: sld_varlist.append( (alphase,4,4.1,alpha) )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = 'OM_rain.in'
sld_varlist = [ ('g2',1) ] 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = 'kinspc.in'
# sld_varlist = [ ('g2',0) ] 
sld_varlist = [ ('g2',g2_turnoverTime) ] 
if include_Al:  sld_varlist.append(  (alphase,0)  )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = '2ndslds.in'
srcfile = os.path.join(datadir,'2ndslds_def.in')
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,srcfile = srcfile
    )





error = 1e4
tol = 1e-4
tol = 1e-3
# tol = 1e-2

cnt = 1
res_list = []

# reading the iterations already done if any
if prev_iter_exist:
    print('previous iteration exists',iter_prev.shape)
    for i in range(iter_prev.shape[0]):
        res_list.append([ int(iter_prev[i,j]) if j == 0 else iter_prev[i,j] for j in range(iter_prev.shape[1])  ])
    
    # cnt,phint_field,phint,omint,acint,error,targetpH,targetOM,acidsat,ca,omrain_field,np.log10(kh)
    cnt = int(iter_prev[-1,0])
    ca  = iter_prev[-1,-3]
    omrain_field = iter_prev[-1,-2]
    kh  = 10.**iter_prev[-1,-1]




# --- define function for a given set 
def run_a_single_set(ca,logkh,omrain_field, rain_list_field):

    # =========== field sim =========== 
    # 4 boundary conditions 
    # (1) ca boundary in M 
    if liming:
        CC_rain = ca
        CC_rain = CC_rain/omrain_field
        rain_list_field = []
    else:
        rain_list_field = rain_list_field # [('ca',ca)]
    
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_field
        ,pr_list = pr_list_field
        ,rain_list=rain_list_field
        ,atm_list=atm_list_field
        )

    # (2) CEC  
    filename = 'cec.in'
    logkhna = logkh
    logkhk  = logkh - 1.1
    logkhca = (logkh - 0.665)*2.
    logkhmg = (logkh - 0.507)*2.
    logkhal = (logkh - 0.41)*3.
    sld_varlist = [('inrt',cec, logkhna, logkhk, logkhca, logkhmg, logkhal, alpha) ,('g2',cec, logkhna, logkhk, logkhca, logkhmg, logkhal, alpha) ] 
    if include_Al:   sld_varlist.append( (alphase,cec,logkhna, logkhk, logkhca, logkhmg, logkhal,alpha) )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # (3) Framework to reflect OM rain
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_field
        ,ztot=ztot_field
        ,nz=nz
        ,ttot=ttot_field
        ,temp=temp_field
        ,fdust=fdust_field
        ,fdust2=fdust2
        ,taudust=taudust_field
        ,omrain=omrain_field
        ,zom=zom
        ,poro=poro_field
        ,moistsrf=moistsrf_field
        ,zwater=zwater
        ,zdust=zdust_field
        ,w=w_field
        ,q=q_field
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_field
        )

    # (4) OM_rain.in to reflect N rain 
    # N_rain = 8.406375  # gN/m2/yr
    N_rain = nitrif  # gN/m2/yr
    N_rain = N_rain/14  # mol N/m2/yr
    N_rain = N_rain*80  # g NH4NO3/m2/yr
    N_rain = N_rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
    N_rain = N_rain/(omrain_field)  # normalize against OC rain 
    
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    if include_N: sld_varlist.append( ('amnt',N_rain) ) 
    if liming: sld_varlist.append( (limesp,CC_rain) )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    
    # >>> run 
    os.system(outdir+runname_field+where+exename)
    
    # >>> 6 data retrieval after field run  
    
    # (1) get porewater pH 
    phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
    
    # (2) get acidity in %
    # dacint = get_int_prof.get_ac_int_site(outdir,runname,dep_sample)
    acint = get_int_prof.get_ac_int_site_v2(outdir,runname_field,dep_sample)
    
    # (3) get bulk density
    # dense = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    
    # (4) get field solid wt%
    sps = sld_list_field # ['g2','inrt']
    if include_Al:   sps.append( alphase ) 
    if liming:   sps.append( limesp ) 
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    # (5) get SOM wt%
    omint = sldwt_list[sps.index('g2')]
    
    # (6) get exchangeable solute concs. 
    aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
    
    # (7) get DIC conc. (added 3.22.2023)
    dic,dep = get_int_prof.get_ave_DIC_save(outdir,runname_field,dep_sample)
    
    print(aqsps,btmconcs,dep)
    
    # =========== lab sim =========== 
    # 2 processes before setting boundary conditions 
    # (1) getting porosity in the laboratory 
    poro_lab = water_frac/(1./dense_lab+water_frac)
    
    # (2) getting material to be dissolved in the laboratory 
    
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
        fdust = ztot_lab*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
    
        if liming and sp == 'ca':  # (added 6.15.2023)
            fdust += ztot_lab*(1-poro_lab)*dense_lab*1e6*sldwt_list[sps.index(limesp)]/100./limemwt * oxide_mass_list[iox]/oxide_stch_list[iox] # (added 6.15.2023)
                
        fdust_list.append(fdust)
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
        if liming and limesp=='cc':  # (added 6.15.2023)
            dic += dense_lab*1e6*sldwt_list[sps.index(limesp)]/100./limemwt # (added 6.15.2023)
            
        fdust_dic = ztot_lab*(1-poro_lab)*dic* 30.      # (added 3.22.2023)
            
        fdust_list.append(fdust_dic)                    # (added 3.22.2023)
        fdust_nm_list.append('g1')                      # (added 3.22.2023)

    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]
    
    # 4 boundary conditions 
    # (1) frame.in to reflect the reference dust flux (fdust_lab)
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_lab
        ,ztot=ztot_lab
        ,nz=nz
        ,ttot=ttot_lab
        ,temp=temp_lab
        ,fdust=fdust_lab
        ,fdust2=fdust2
        ,taudust=taudust_lab
        ,omrain=omrain_lab
        ,zom=zom
        ,poro=poro_lab
        ,moistsrf=moistsrf_lab
        ,zwater=zwater
        ,zdust=zdust_lab
        ,w=w_lab
        ,q=q_lab
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_lab
        )
    
    # (2) solid phase concs. 
    if "inrt" in sldwt_list:
        pr_list_lab = [('inrt',sldwt_list[sps.index('inrt')]/100.),('g2',sldwt_list[sps.index('g2')]/100.)]
    else:
        pr_list_lab = [('g2',sldwt_list[sps.index('g2')]/100.)]
    if include_Al:  pr_list_lab.append(  (alphase,sldwt_list[sps.index(alphase)]/100.) )
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )

    # (3) CEC 
    filename = 'cec.in'
    sld_varlist = [('inrt',cec, logkhna, logkhk, logkhca, logkhmg, logkhal,alpha),('g2',cec,logkhna, logkhk, logkhca, logkhmg, logkhal,alpha)  ] 
    if include_Al:   sld_varlist.append( (alphase,cec,logkhna, logkhk, logkhca, logkhmg, logkhal,alpha) )
    if include_DIC:  sld_varlist.append( ('g1',0, logkhna, logkhk, logkhca, logkhmg, logkhal,alpha) )    # added 3.23.2023
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # (4) dust composition
    filename = 'dust.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # >>> run 
    os.system(outdir+runname_lab+where+exename)

    # 1 data retrieval 
    # (1) lab pH
    phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
    
    return phint_field,phint,omint,acint






# --- RUN A SINGLE RUN
while (error > tol):

    ymx = np.zeros(3)
    emx = np.zeros(3)
    amx = np.zeros((3,3),dtype=np.float64)
    
    facts = [1e-4]*3
    # facts = [1e-8,1e-2,1e-4]
    # facts = [1e-2]*3

    # command = ' {:.8f} {:.8f} {:.8f} {}'.format(cec,np.log10(kh),ca,runname)  
    
    # =========== sim 1st =========== 
    
    phint_field,phint,omint,acint = run_a_single_set(ca,np.log10(kh),omrain_field, rain_list_field)

    print(phint_field,phint,acint)

    # # =========== sim 2nd =========== 
    # dca = facts[0]
    # dca = ca*facts[0]
    
    # dphint_field,dphint,domint,dacint = run_a_single_set(ca+dca,np.log10(kh),omrain_field)
    
    # if phnorm_pw:       dphint_dca = (dphint_field-phint_field)/dca * ca
    # if not phnorm_pw:   dphint_dca = (dphint-phint)/dca * ca
    # # if phnorm_pw:       dphint_dca = (10**-dphint_field-10**-phint_field)/dca * ca
    # # if not phnorm_pw:   dphint_dca = (10**-dphint-10**-phint)/dca * ca
    # dacint_dca = (dacint-acint)/dca * ca
    # domint_dca = (domint-omint)/dca * ca

    # print(dphint_dca,dacint_dca,domint_dca)

    # # =========== sim 3rd =========== 
    
    # dkh = facts[1]
    # dkh = kh*facts[1]
    
    # dphint_field,dphint,domint,dacint = run_a_single_set(ca,np.log10(kh+dkh),omrain_field)

    # # domint_dlogkh = (domint-omint)/dlogkh
    # # dacint_dlogkh = (dacint-acint)/dlogkh
    # # if phnorm_pw:       dphint_dlogkh = (dphint_field-phint_field)/dlogkh
    # # if not phnorm_pw:   dphint_dlogkh = (dphint-phint)/dlogkh 
    # # if phnorm_pw:       dphint_dlogkh = (10**-dphint_field-10**-phint_field)/dkh * kh
    # # if not phnorm_pw:   dphint_dlogkh = (10**-dphint-10**-phint)/dkh * kh
    # if phnorm_pw:       dphint_dlogkh = (dphint_field-phint_field)/dkh*kh
    # if not phnorm_pw:   dphint_dlogkh = (dphint-phint)/dkh*kh 
    # dacint_dlogkh = (dacint-acint)/dkh * kh
    # domint_dlogkh = (domint-omint)/dkh * kh

    # print(dphint_dlogkh,dacint_dlogkh,domint_dlogkh)
    

    # # =========== sim 4th =========== 
    # domrain = facts[2]
    # domrain = omrain_field*facts[2]
    
    # dphint_field,dphint,domint,dacint = run_a_single_set(ca,np.log10(kh),omrain_field+domrain)

    # if phnorm_pw:       dphint_domrain = (dphint_field-phint_field)/domrain * omrain_field
    # if not phnorm_pw:   dphint_domrain = (dphint-phint)/domrain * omrain_field
    # # if phnorm_pw:       dphint_domrain = (10**-dphint_field-10**-phint_field)/domrain * omrain_field
    # # if not phnorm_pw:   dphint_domrain = (10**-dphint-10**-phint)/domrain * omrain_field
    # dacint_domrain = (dacint-acint)/domrain * omrain_field
    # domint_domrain = (domint-omint)/domrain * omrain_field

    # print(dphint_domrain,dacint_domrain,domint_domrain)
    
    # # ===========  Newton method ==================
    # # filling Jacobian matrix 

    if phnorm_pw:       ymx[0] = phint_field - targetpH
    if not phnorm_pw:   ymx[0] = phint - targetpH
    # if phnorm_pw:       ymx[0] = 10**-phint_field - 10**-targetpH
    # if not phnorm_pw:   ymx[0] = 10**-phint - 10**-targetpH
    ymx[1] = (acint - acidsat) 
    ymx[2] = omint - targetOM 

    ymx = -ymx

    print("----------------------------------------- ")
    print("")
    print("ph target: " + str(targetpH) + " | ph offset: " + str(ymx[0]))
    print("")
    print("----------------------------------------- ")
    
    time.sleep(5)
    
    res_list.append([cnt,phint_field,phint,omint,acint,error,targetpH,targetOM,acidsat,ca,omrain_field,np.log10(kh), g2_turnoverTime])
    cnt += 1
    
    if cnt > iter_max: break
    
    name_list = [
        'iter.'
        ,'porewater_pH[-]'
        ,'soil_pHw[-]'
        ,'OM[wt%]'
        ,'exchangeable_acidity[%CEC]'
        ,'error'
        ,'target_pH[-]'
        ,'target_OM[wt%]'
        ,'target_exchangeable_acidity[%CEC]'
        ,'Ca[M]' if not liming else 'liming[g:{}/m2/yr]'.format(limesp.upper())
        ,'OM_rain[gC/m2/yr]'
        ,'log10(KH/Na)'
        ]
    for runname in [runname_field,runname_lab]:
        # np.savetxt(outdir + runname + where + 'iteration_tmp.res',np.array(res_list))
        dst = outdir + runname + where + 'iteration_tmp.res'

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

    error = 1e-99

    
    if use_local_storage:
        dst2 = './Newton/'+ simid + '_iteration_tmp.res'
        shutil.copyfile(dst,dst2)





# --- PRINT SOME DIAGNOSTICS
# 
phint_field,phint,omint,acint = run_a_single_set(ca,np.log10(kh),omrain_field, rain_list_field)
print(phint_field,phint,omint,acint)

print("----------------------------------------- ")
print("")
print("ph target: " + str(targetpH) + " | ph offset: " + str(ymx[0]))
print("")
print("----------------------------------------- ")

# name_list = [
    # 'iter.'
    # ,'porewater_pH[-]'
    # ,'soil_pHw[-]'
    # ,'error'
    # ]
name_list = [
    'iter.',
    'porewater_pH[-]',
    'soil_pHw[-]',
    'OM[wt%]',
    'exchangeable_acidity[%CEC]',
    'error',
    'target_pH[-]',
    'target_OM[wt%]',
    'target_exchangeable_acidity[%CEC]',
    'Ca[M]' if not liming else 'liming[g{:}/m2/yr]'.format(limesp.upper()),
    'OM_rain[gC/m2/yr]',
    'log10(KH/Na)',
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


# ... move directory to AWS
if aws_save == "move":
    for runname in [runname_field,runname_lab]:
        src = os.path.join(outdir, runname)
        dst_aws = os.path.join(aws_bucket)
        # check if the destination directory exists and remove it if it does
        # if os.path.exists(dst_aws):
        #     shutil.rmtree(dst_aws)
        # move the dir
        # shutil.move(src, dst_aws)
        cmd_activate = 'conda run -n cdr-scepter1p0_env'
        cmd_run = 's5cmd mv ' + src + ' ' + dst_aws
        result1 = subprocess.run(cmd_activate, shell=True, check=True)
        result2 = subprocess.run(cmd_run, shell=True, check=True)
        #print(result)

elif aws_save == "copy":
    for runname in [runname_field,runname_lab]:
        src = os.path.join(outdir, runname)
        dst_aws = os.path.join(aws_bucket)
        # check if the destination directory exists and remove it if it does
        # if os.path.exists(dst_aws):
            # shutil.rmtree(dst_aws)
        # copy the dir
        # shutil.copytree(src, dst_aws)
        cmd_activate = 'conda run -n cdr-scepter1p0_env'
        cmd_run = 's5cmd cp ' + src + ' ' + dst_aws
        result1 = subprocess.run(cmd_activate, shell=True, check=True)
        result2 = subprocess.run(cmd_run, shell=True, check=True)
        #print(result)
        
    