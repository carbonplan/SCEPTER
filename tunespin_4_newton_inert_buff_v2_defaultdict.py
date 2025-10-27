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

# --- read in helper functions from ew_workflows
# import module
from ew_workflows import scepter_helperFxns as shf
# ---


# -------------------------------------------------------------
# -------------------------------------------------------------
# --- set default and system args
sys_args = shf.parse_arguments(sys.argv)   # parse system args
import_dict = sys_args['default_dict'] # set the dictionary to use from system args
def_args = getattr(defaults.dict_singlerun, import_dict)  # get dict attribute

# set global variables
# (TK: currently, combined_dict doesn't get saved in this version)
combined_dict = shf.set_vars(def_args, sys_args)  # default unless defined in sys_args
# (set as kwargs CHANGE WHEN THIS BECOMES A FUNCTION)
kwargs = combined_dict.copy()
# add to globals (CHANGE/REMOVE WHEN THIS BECOMES A FUNCTION
for key, value in combined_dict.items():
    globals()[key] = value
# (vars are saved lower after we make the directory)

# --- update optional inputs
ttot_field = 10000 if kwargs.get('ttot_field') is None else kwargs.get('ttot_field')
ztot_field = 0.5 if kwargs.get('ztot_field') is None else kwargs.get('ztot_field')
zom = 0.25 if kwargs.get('zom') is None else kwargs.get('zom')
if not make_initial_guess:
    omrain_field = 900 if kwargs.get('omrain_field') is None else kwargs.get('omrain_field')
zwater = 10000 if kwargs.get('zwater') is None else kwargs.get('zwater')
w_scheme_field = 1 if kwargs.get('w_scheme_field') is None else kwargs.get('w_scheme_field')
mix_scheme_field = 1 if kwargs.get('mix_scheme_field') is None else kwargs.get('mix_scheme_field')
poro_iter_field = 'false' if kwargs.get('poro_iter_field') is None else kwargs.get('poro_iter_field')
poro_evol = 'false' if kwargs.get('poro_evol') is None else kwargs.get('poro_evol')
sldmin_lim = 'false' if kwargs.get('sldmin_lim') is None else kwargs.get('sldmin_lim')
psd_bulk_field = 'false' if kwargs.get('psd_bulk_field') is None else kwargs.get('psd_bulk_field')
psd_full_field = 'false' if kwargs.get('psd_full_field') is None else kwargs.get('psd_full_field')
sa_evol_1 = 'true' if kwargs.get('sa_evol_1') is None else kwargs.get('sa_evol_1')
sa_evol_2 = 'false' if kwargs.get('sa_evol_2') is None else kwargs.get('sa_evol_2')
display = 'true' if kwargs.get('display') is None else kwargs.get('display')
disp_lim = 'true' if kwargs.get('disp_lim') is None else kwargs.get('disp_lim')
close_aq_field = 'false' if kwargs.get('close_aq_field') is None else kwargs.get('close_aq_field')
season = 'false' if kwargs.get('season') is None else kwargs.get('season') 
nz = 30 if kwargs.get('nz') is None else int(kwargs.get('nz'))

water_frac = water_frac_tunespin

targetpH = tph
acidsat = tec
targetOM = tsom
targetlogpco2 = tsoilco2
poro_input = poro
moist_input = soilmoisture
sat_input = moist_input/poro
kh = 10**(3.5 + alpha/2.)

# tau_g2 = 8. # yr, turn over time for OM
tau_g2 = 1. # yr, turn over time for OM
tau_g2 = 1. / ( 1./8.*3.**( (mat - 15.)/10.  ) )
# print(tau_g2)

datadir = os.path.join(modeldir, 'data/')

runname_guess = initguess

if runname_guess=='none':
    make_initial_guess = False
else:
    make_initial_guess = True

if make_initial_guess:
    outdir_guess = outdir
    print("my initial guess ----")
    try:
        filename = 'iteration.res'  # [tk] deleted leading "/" because it doesn't play nicely with "os.path.join()"
        try:
            filename = 'iteration.res'  # [tk] deleted leading "/" because it doesn't play nicely with "os.path.join()"
            if 'output' in outdir_guess:
                data_tmp = np.loadtxt(os.path.join(outdir_guess, runname_guess, filename),skiprows=1)
            if 'archive' in outdir_guess:
                tar_tmp = tarfile.open(os.path.join(outdir_guess, runname_guess) + '.tar.gz',"r")
                xfile = tar_tmp.extractfile('./'+runname_guess + filename,skiprows=1)
                data_tmp = np.loadtxt(xfile)
        except:
            filename = 'iteration_tmp.res' # 'DA_iteration_tmp.res'
            if 'output' in outdir_guess:
                data_tmp = np.loadtxt(os.path.join(outdir_guess, runname_guess, filename),skiprows=1)
            if 'archive' in outdir_guess:
                tar_tmp = tarfile.open(outdir_guess + runname_guess + '.tar.gz',"r")
                xfile = tar_tmp.extractfile('./'+runname_guess + filename,skiprows=1)
                data_tmp = np.loadtxt(xfile)
        # except:
            # filename = '/DA_iteration_tmp.res'
            # if 'output' in outdir_guess:
                # data_tmp = np.loadtxt(outdir_guess + runname_guess + filename)
            # if 'archive' in outdir_guess:
                # tar_tmp = tarfile.open(outdir_guess + runname_guess + '.tar.gz',"r")
                # xfile = tar_tmp.extractfile('./'+runname_guess + filename)
                # data_tmp = np.loadtxt(xfile)
        if len(data_tmp.shape)==2:
            j_col = 1
        else:
            j_col = 0

        if data_tmp.shape[j_col]==12:
            i_error = 5
            i_log10kh = -1
            i_omrain_field = -2
            i_ca = -3
        else:
            i_error = 6
            i_log10kh = -2
            i_omrain_field = -3
            i_ca = -4
            i_tau = -1

        try:
            i_errormin = np.argmin( data_tmp[:,i_error] )
            error_min = data_tmp[i_errormin,i_error]*100.
            log10kh = (data_tmp[i_errormin,i_log10kh])
            omrain_field = (data_tmp[i_errormin,i_omrain_field])
            ca = (data_tmp[i_errormin,i_ca])
        except:
            error_min = data_tmp[i_error]*100.
            log10kh = (data_tmp[i_log10kh])
            omrain_field = (data_tmp[i_omrain_field])
            ca = (data_tmp[i_ca])
        
        filename = 'cec.in'
        sld_data_list = get_inputs.get_input_sld_properties(outdir_guess,runname_guess,filename)
        alpha_ref = sld_data_list[0][-1]
        
        kh = 10**( (alpha - alpha_ref)/2. + log10kh )
        
        print( 'initial guess from {:}'.format(runname_guess) )
        print( np.log10(ca), np.log10(kh), omrain_field )
        
        if data_tmp.shape[j_col]>12:
            try:
                i_errormin = np.argmin( data_tmp[:,i_error] )
                tau_g2 = (data_tmp[i_errormin,i_tau])
            except:
                tau_g2 = (data_tmp[i_tau])

            print( 'additional initial guess from {:}'.format(runname_guess) )
            print( tau_g2 )
    except:
        print( 'initial guess from {:} not used as it did not converge'.format(runname_guess) )
        make_initial_guess = False
        if stop_unsuccessful: exit()
    
    # exit()
import time
time.sleep(10)

if liming: ca = 10


# save file denoting the variables used
outdir = outdir
if use_local_storage:  outdir = os.environ['TMPDIR'] + '/scepter_output/'
simid = spinname
runname_field   = simid+'_spintuneup4_field'
runname_lab     = simid+'_spintuneup4_lab'



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
# os.system('make --file=makefile_test')
prev_iter_exist = False
if not os.path.exists( outdir + runname_field) : 
    os.system('mkdir -p ' + outdir + runname_field)
else:
    if os.path.exists( outdir + runname_field + where + 'iteration_tmp.res'):
        try:
            iter_prev = np.loadtxt(outdir + runname_field + where + 'iteration_tmp.res',skiprows=1)
            shutil.copy(outdir + runname_field + where + 'iteration_tmp.res'
                ,outdir + runname_field + where + 'iteration_tmp_SAVE_'+datestr+'.res')
        except:
            print('something went wrong while saving previous iterations')
        # prev_iter_exist = True
if not os.path.exists( outdir + runname_lab) : os.system('mkdir -p ' + outdir + runname_lab)

# save the vars used ---
fn_dict_save = os.path.join(outdir, runname_field, "vars.res")
shf.save_dict_to_text_file(combined_dict, fn_dict_save, delimiter='\t')
# ----------------------

os.system('cp ' + exename_src + to + outdir + runname_field + where + exename)
os.system('cp ' + exename_src + to + outdir + runname_lab + where + exename)
# change back
os.chdir(mycwd)

ztot=0.5
ztot_lab=0.05
ttot_lab=100
temp_field=mat
temp_lab=25
fdust_field=0
fdust_lab=0
fdust2=0
taudust_field=0
taudust_lab=0.01
omrain_lab=0
# poro_field=0.5
poro_field=poro_input
poro_lab=0.928391508
# moistsrf_field=0.5
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


restart ='false'
if include_roughness_sa == True:
    rough_field      ='true'
else:
    rough_field      ='false'
act_ON ='false'
if activity_on : act_ON ='true'
dt_fix='false'
if cec_adsorption_on == True:
    cec_on="true"
else:
    cec_on="false"
dz_fix='true'

w_scheme_lab=0
mix_scheme_lab=0 
poro_iter_lab='true' 
rough_lab      ='false'
close_aq_lab='true'
psd_bulk_lab='false'
psd_full_lab ='false'
    
make_inputs.get_input_switches(
    outdir=outdir
    ,runname=runname_field
    ,w_scheme=w_scheme_field
    ,mix_scheme=mix_scheme_field
    ,poro_iter=poro_iter_field
    ,sldmin_lim=sldmin_lim 
    ,display=display
    # ,report=report
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

sld_list_field=['inrt','g2']
if liming:
    sld_list_field.append(limesp)
aq_list_field = ['ca','k','mg','na']
if include_N: 
    sld_list_field.append('amnt')
    aq_list_field.append('no3')
if include_Al: 
    sld_list_field.append(alphase)
    aq_list_field.append('al')
if add_secondary:                # [tk added to track secondary precipitates]
    sld_list_field += sld_track

gas_list_field = ['pco2']
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

# bring in parent rock file [tk added]
if (spinup_parentrock_file != None) & ("json" in spinup_parentrock_file):
    pr_file = os.path.join(modeldir, 'data', spinup_parentrock_file)
    with open(pr_file, 'r') as f:
        pr_list_field = json.load(f)
else:
    pr_list_field = [('inrt',1.0)]
    
rain_list_field = [('ca',5.0e-6)]
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
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,srcfile = srcfile
    )
    
filename = 'kinspc.in'
sld_varlist = [ ('g2',tau_g2) ] 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
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
    # ,report=report
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
sld_varlist = [ ('g2',0) ] 
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
    

def run_a_single_set(ca,logkh,omrain_field,tau_g2):

    # =========== field sim =========== 
    # 4 boundary conditions 
    # (1) ca boundary in M 
    if liming:
        CC_rain = ca
        CC_rain = CC_rain/omrain_field
        rain_list_field = []
    else:
        rain_list_field = [('ca',ca)]
    
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
    
    # (5) OM decay const | added 06-16-2023
    filename = 'kinspc.in'
    sld_varlist = [ ('g2',tau_g2) ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    
    # >>> run 
    # os.system(outdir+runname_field+where+exename)
    os.system("chmod +x " + os.path.join(outdir,runname_field,'scepter'))  # grant permissions
    proc = subprocess.Popen([outdir+runname_field+where+exename])

    my_timeout = tunespin_timeout # 120*60  # kill the run if not complete in this many seconds
    
    flag_error_field = False
    flag_error_lab = False
    
    try:
        proc.wait(my_timeout)
        print('run finished within {:f} min'.format(int(my_timeout/60.)))
    except subprocess.TimeoutExpired:
        proc.kill()
        flag_error_field = True
        phint_field=np.nan;phint=np.nan;omint=np.nan;acint=np.nan;logpco2=np.nan
        print('run UNfinished within {:f} min'.format(int(my_timeout/60.)))
        
        if stop_unsuccessful: exit()
        
        return phint_field,phint,omint,acint,logpco2,flag_error_field,flag_error_lab
    
    # >>> 6 data retrieval after field run  
    
    flag_error_field = False
    flag_error_lab = False
    
    try:
        # (1) get porewater pH 
        # phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample,20)
        phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample) # [tykukla] rm 20
        
        # (2) get acidity in %
        # dacint = get_int_prof.get_ac_int_site(outdir,runname,dep_sample)
        # acint = get_int_prof.get_ac_int_site_v2(outdir,runname_field,dep_sample,20)
        acint = get_int_prof.get_ac_int_site_v2(outdir,runname_field,dep_sample) # [tykukla] rm 20
        
        # (3) get bulk density
        # dense = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
        # dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample,20)
        dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample) # [tykukla] rm 20
        
        # (4) get field solid wt%
        sps = ['g2','inrt']
        if include_Al:   sps.append( alphase ) 
        if liming:   sps.append( limesp ) 
        sldwt_list = []
        for sp in sps:
            # sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp],20)
            sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp]) # [tykukla] rm 20
            sldwt_list.append(sldwt)
        exchanger = sum(sldwt_list)
        
        # (5) get SOM wt%
        omint = sldwt_list[sps.index('g2')]
        
        # (6) get exchangeable solute concs. 
        # aqsps,btmconcs,dep,time = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample,20)  # returning mol/ solid m3 depth averaged value 
        aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample) # [tykukla] rm 20 and delete time from results (not output in this version of .get_totsave_..) ;  # returning mol/ solid m3 depth averaged value 
        
        # (7) get DIC conc. (added 3.22.2023)
        # dic,dep = get_int_prof.get_ave_DIC_save(outdir,runname_field,dep_sample,20)
        dic,dep = get_int_prof.get_ave_DIC_save(outdir,runname_field,dep_sample) # [tykukla] rm 20
        
        # (8) get soil pCO2 (added 6.16.2023)
        # spco2 = get_int_prof.get_gas_int_site(outdir,runname_field,dep_sample,'pco2',20)
        spco2 = get_int_prof.get_gas_int_site(outdir,runname_field,dep_sample,'pco2') # [tykukla] rm 20
        logpco2 = np.log10(spco2)
        
        print(aqsps,btmconcs,dep)
        print('!!! SUCCESSFUL data retrieval; field run must have been SUCCESSFUL !!!')
    except:
        print('*** UNSUCCESSFUL data retrieval; field run must have been UNSUCCESSFUL ***')
        flag_error_field = True
        phint_field=np.nan;phint=np.nan;omint=np.nan;acint=np.nan;logpco2=np.nan
        
        if stop_unsuccessful: exit()
        
        return phint_field,phint,omint,acint,logpco2,flag_error_field,flag_error_lab
    
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
    pr_list_lab = [('inrt',sldwt_list[sps.index('inrt')]/100.),('g2',sldwt_list[sps.index('g2')]/100.)]
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
        
    # before running previous results are deleted
    if os.path.exists(outdir+runname_lab+where+'flx'): shutil.rmtree(outdir+runname_lab+where+'flx')
    if os.path.exists(outdir+runname_lab+where+'prof'): shutil.rmtree(outdir+runname_lab+where+'prof')
    
    # >>> run 
    # os.system(outdir+runname_lab+where+exename)
    os.system("chmod +x " + os.path.join(outdir,runname_lab,'scepter'))  # grant permissions
    proc = subprocess.Popen([outdir+runname_lab+where+exename])

    my_timeout = tunespin_timeout # 60*20
    
    flag_error_field = False
    flag_error_lab = False
    
    try:
        proc.wait(my_timeout)
        print('run finished within {:f} min'.format(int(my_timeout/60.)))
    except subprocess.TimeoutExpired:
        proc.kill()
        flag_error_lab = True
        phint_field=np.nan;phint=np.nan;omint=np.nan;acint=np.nan;logpco2=np.nan
        print('run UNfinished within {:f} min'.format(int(my_timeout/60.)))
        if stop_unsuccessful: exit()
        return phint_field,phint,omint,acint,logpco2,flag_error_field,flag_error_lab

    
    # 1 data retrieval 

    flag_error_field = False
    flag_error_lab = False
    try:
        # (1) lab pH
        # phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample,20)
        phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)  # [tykukla] delete 20
        print('!!! SUCCESSFUL data retrieval; lab run must have been SUCCESSFUL !!!')
    
    except:
        print('*** UNSUCCESSFUL data retrieval; lab run must have been UNSUCCESSFUL ***')
        flag_error_lab = True
        phint_field=np.nan;phint=np.nan;omint=np.nan;acint=np.nan;logpco2=np.nan
        if stop_unsuccessful: exit()
        return phint_field,phint,omint,acint,logpco2,flag_error_field,flag_error_lab
    
    return phint_field,phint,omint,acint,logpco2,flag_error_field,flag_error_lab

while (error > tol):

    nmx = 4

    ymx = np.zeros(nmx)
    emx = np.zeros(nmx)
    amx = np.zeros((nmx,nmx),dtype=np.float64)
    
    facts = [1e-4]*nmx
    # facts = [1e-8,1e-2,1e-4]
    # facts = [1e-2]*3

    # command = ' {:.8f} {:.8f} {:.8f} {}'.format(cec,np.log10(kh),ca,runname)  
    
    # =========== sim 1st =========== 
    
    phint_field,phint,omint,acint,logpco2,flag_error_field,flag_error_lab = run_a_single_set(ca,np.log10(kh),omrain_field,tau_g2)
    
    if flag_error_field or flag_error_lab: # assume error is caused by too much OM rain
        print(flag_error_field,flag_error_lab,omrain_field)
        # omrain_field = omrain_field/1.5
        omrain_field = omrain_field/10.
        continue

    print('after 1st run',phint_field,phint,acint,logpco2)
    time.sleep(5)

    # =========== sim 2nd =========== 
    dca = facts[0]
    dca = ca*facts[0]
    
    dphint_field,dphint,domint,dacint,dlogpco2,flag_error_field,flag_error_lab = run_a_single_set(ca+dca,np.log10(kh),omrain_field,tau_g2)
    print('after 2nd run',dphint_field,dphint,domint,dacint,dlogpco2)
    time.sleep(5)
    
    if phnorm_pw:       dphint_dca = (dphint_field-phint_field)/dca * ca
    if not phnorm_pw:   dphint_dca = (dphint-phint)/dca * ca
    # if phnorm_pw:       dphint_dca = (10**-dphint_field-10**-phint_field)/dca * ca
    # if not phnorm_pw:   dphint_dca = (10**-dphint-10**-phint)/dca * ca
    dacint_dca = (dacint-acint)/dca * ca
    domint_dca = (domint-omint)/dca * ca
    dlogpco2_dca = (dlogpco2-logpco2)/dca * ca


    # =========== sim 3rd =========== 
    
    dkh = facts[1]
    dkh = kh*facts[1]
    
    dphint_field,dphint,domint,dacint,dlogpco2,flag_error_field,flag_error_lab = run_a_single_set(ca,np.log10(kh+dkh),omrain_field,tau_g2)

    # domint_dlogkh = (domint-omint)/dlogkh
    # dacint_dlogkh = (dacint-acint)/dlogkh
    # if phnorm_pw:       dphint_dlogkh = (dphint_field-phint_field)/dlogkh
    # if not phnorm_pw:   dphint_dlogkh = (dphint-phint)/dlogkh 
    # if phnorm_pw:       dphint_dlogkh = (10**-dphint_field-10**-phint_field)/dkh * kh
    # if not phnorm_pw:   dphint_dlogkh = (10**-dphint-10**-phint)/dkh * kh
    if phnorm_pw:       dphint_dlogkh = (dphint_field-phint_field)/dkh*kh
    if not phnorm_pw:   dphint_dlogkh = (dphint-phint)/dkh*kh 
    dacint_dlogkh = (dacint-acint)/dkh * kh
    domint_dlogkh = (domint-omint)/dkh * kh
    dlogpco2_dlogkh = (dlogpco2-logpco2)/dkh * kh

    print('after 3rd run',dphint_dlogkh,dacint_dlogkh,domint_dlogkh,dlogpco2_dlogkh)
    time.sleep(5)
    

    # =========== sim 4th =========== 
    domrain = facts[2]
    domrain = omrain_field*facts[2]
    
    dphint_field,dphint,domint,dacint,dlogpco2,flag_error_field,flag_error_lab = run_a_single_set(ca,np.log10(kh),omrain_field+domrain,tau_g2)

    if phnorm_pw:       dphint_domrain = (dphint_field-phint_field)/domrain * omrain_field
    if not phnorm_pw:   dphint_domrain = (dphint-phint)/domrain * omrain_field
    # if phnorm_pw:       dphint_domrain = (10**-dphint_field-10**-phint_field)/domrain * omrain_field
    # if not phnorm_pw:   dphint_domrain = (10**-dphint-10**-phint)/domrain * omrain_field
    dacint_domrain = (dacint-acint)/domrain * omrain_field
    domint_domrain = (domint-omint)/domrain * omrain_field
    dlogpco2_domrain = (dlogpco2-logpco2)/domrain * omrain_field

    print('after 4th run',dphint_domrain,dacint_domrain,domint_domrain,dlogpco2_domrain)
    time.sleep(5)
    

    # =========== sim 5th =========== 
    dtaug2 = facts[3]
    dtaug2 = tau_g2*facts[3]
    
    dphint_field,dphint,domint,dacint,dlogpco2,flag_error_field,flag_error_lab = run_a_single_set(ca,np.log10(kh),omrain_field,tau_g2+dtaug2)

    if phnorm_pw:       dphint_dtaug2 = (dphint_field-phint_field)/dtaug2 * tau_g2
    if not phnorm_pw:   dphint_dtaug2 = (dphint-phint)/dtaug2 * tau_g2
    # if phnorm_pw:       dphint_domrain = (10**-dphint_field-10**-phint_field)/domrain * omrain_field
    # if not phnorm_pw:   dphint_domrain = (10**-dphint-10**-phint)/domrain * omrain_field
    dacint_dtaug2 = (dacint-acint)/dtaug2 * tau_g2
    domint_dtaug2 = (domint-omint)/dtaug2 * tau_g2
    dlogpco2_dtaug2 = (dlogpco2-logpco2)/dtaug2 * tau_g2

    print('after 5th run',dphint_dtaug2,dacint_dtaug2,domint_dtaug2,dlogpco2_dtaug2)
    time.sleep(5)
    
    # ===========  Newton method ==================
    # filling Jacobian matrix 

    if phnorm_pw:       ymx[0] = phint_field - targetpH
    if not phnorm_pw:   ymx[0] = phint - targetpH
    # if phnorm_pw:       ymx[0] = 10**-phint_field - 10**-targetpH
    # if not phnorm_pw:   ymx[0] = 10**-phint - 10**-targetpH
    ymx[1] = (acint - acidsat) 
    ymx[2] = omint - targetOM 
    ymx[3] = logpco2 - targetlogpco2 

    amx[0,0] = dphint_dca
    amx[0,1] = dphint_dlogkh
    amx[0,2] = dphint_domrain
    amx[0,3] = dphint_dtaug2
    
    amx[1,0] = dacint_dca
    amx[1,1] = dacint_dlogkh
    amx[1,2] = dacint_domrain
    amx[1,3] = dacint_dtaug2
    
    amx[2,0] = domint_dca
    amx[2,1] = domint_dlogkh
    amx[2,2] = domint_domrain
    amx[2,3] = domint_dtaug2
    
    amx[3,0] = dlogpco2_dca
    amx[3,1] = dlogpco2_dlogkh
    amx[3,2] = dlogpco2_domrain
    amx[3,3] = dlogpco2_dtaug2

    ymx = -ymx

    dx = np.linalg.solve(amx, ymx)

    print(ca,np.log10(kh),omrain_field,tau_g2)
    print(ca*np.exp(dx[0]),np.log10(kh*np.exp(dx[1])),omrain_field*np.exp(dx[2]),tau_g2*np.exp(dx[3]) )

    if np.isnan(ymx).any():
        print('nan detected in solution')
        error = 1e-99
    else:
        thrs = [5]*nmx
        # thrs = [3]*3
        # thrs = [7]*3
        if dx[0]>thrs[0]:
            ca = 1.5*ca
        elif dx[0]<-1*thrs[0]:
            ca = ca/1.5
        else:
            ca = ca*np.exp(dx[0])
        
        if dx[1] > thrs[1]:
            kh = 1.5*kh
        elif dx[1] < -1*thrs[1]:
            kh = kh/1.5
        else:
            kh = kh*np.exp(dx[1])
        
        if dx[2] > thrs[2]:
            omrain_field = 1.5*omrain_field
        elif dx[2] < -1*thrs[2]:
            omrain_field = omrain_field/1.5
        else:
            omrain_field = omrain_field*np.exp(dx[2])
        
        if dx[3] > thrs[3]:
            tau_g2 = 1.5*tau_g2
        elif dx[3] < -1*thrs[3]:
            tau_g2 = tau_g2/1.5
        else:
            tau_g2 = tau_g2*np.exp(dx[3])
        
        if phnorm_pw:       emx[0] = ymx[0]/targetpH
        if not phnorm_pw:   emx[0] = ymx[0]/targetpH
        # if phnorm_pw:       emx[0] = ymx[0]/10**-targetpH
        # if not phnorm_pw:   emx[0] = ymx[0]/10**-targetpH
        
        if acidsat!=0: 
            emx[1] = ymx[1]/acidsat 
        else:
            if acint!=0: 
                emx[1] = ymx[1]/acint
            else:
                emx[1] = ymx[1]
                # emx[1] = ymx[1]*1e2
                
                
        if targetOM!=0:
            emx[2] = ymx[2]/targetOM 
        else:
            if omrain_field!=0:
                emx[2] = ymx[2]/omrain_field
            else:
                emx[2] = ymx[2]
                
        if targetlogpco2!=0:
            emx[3] = ymx[3]/targetlogpco2 
        else:
            if logpco2!=0:
                emx[3] = ymx[3]/logpco2
            else:
                emx[3] = ymx[3]
                # emx[2] = ymx[2]*1e2

        error = np.max(np.abs(emx))
        
        
        # when log10 kh is below 1, calculation can become difficult
        if kh < 10**1:
            kh = 10**1
            error = 100
        
        if omrain_field > 1e5 and tau_g2>1:
            omrain_field = omrain_field/10.
            
    
    print(' ')
    print(' ')
    print(' ')
    print(' ')
    print('******* error = ',error)
    print(' ')
    print(' ')
    print(' ')
    print(' ')
    
    time.sleep(5)
    
    res_list.append([cnt,phint_field,phint,omint,acint,logpco2,error,targetpH,targetOM,acidsat,targetlogpco2,ca,omrain_field,np.log10(kh),tau_g2])
    cnt += 1

    print(" --------------------------------------- ")
    print(res_list)
    print(" --------------------------------------- ")    
    
    if cnt > iter_max: break
    
    name_list = [
        'n',
        'p.w.pH',
        'soil.pH',
        'OM[wt%]',
        'ex.acd[%]',
        'logpCO2[atm]',
        'error',
        'tgt.pH',
        'tgt.OM[wt%]',
        'tgt.ex.acd[%]',
        'tgt.logpCO2[atm]',
        'Ca[M]' if not liming else 'lime[g:{}/m2/yr]'.format(limesp.upper()),
        'OM.rain[gC/m2/yr]',
        'log10(KH/Na)',
        'tau.OM[yr]',
        ]
    for runname in [runname_field,runname_lab]:
        # np.savetxt(outdir + runname + where + 'iteration_tmp.res',np.array(res_list))
        dst = os.path.join(outdir, runname, 'iteration_tmp.res')

        with open(dst, 'w') as file:
            for item in name_list:
                if name_list.index(item)==0:
                    file.write('{:5}\t'.format(item))
                elif name_list.index(item)==len(name_list)-1:
                    file.write('{:15}\n'.format(item))
                else:
                    file.write('{:15}\t'.format(item))
            for j in range(len(res_list)):
                item_list = res_list[j]
                for i in range(len(item_list)):
                    item = item_list[i]
                    if i==0:
                        file.write('{:5d}\t'.format(item))
                    elif i==len(item_list)-1:
                        file.write('{:.9e}\n'.format(item))
                    else:
                        file.write('{:.9e}\t'.format(item))
                
        print(res_list)
        
    if use_local_storage:
        dst2 = './Newton/'+ simid + '_iteration_tmp.res'
        shutil.copyfile(dst,dst2)


phint_field,phint,omint,acint,logpco2,flag_error_field,flag_error_lab = run_a_single_set(ca,np.log10(kh),omrain_field,tau_g2)
print(phint_field,phint,omint,acint,logpco2)
# name_list = [
    # 'iter.'
    # ,'porewater_pH[-]'
    # ,'soil_pHw[-]'
    # ,'error'
    # ]
    
name_list = [
    'n',
    'p.w.pH',
    'soil.pH',
    'OM[wt%]',
    'ex.acd[%]',
    'logpCO2[atm]',
    'error',
    'tgt.pH',
    'tgt.OM[wt%]',
    'tgt.ex.acd[%]',
    'tgt.logpCO2[atm]',
    'Ca[M]' if not liming else 'lime[g:{}/m2/yr]'.format(limesp.upper()),
    'OM.rain[gC/m2/yr]',
    'log10(KH/Na)',
    'tau.OM[yr]',
    ]
    
for runname in [runname_field,runname_lab]:
    dst = outdir + runname + where + 'iteration.res'

    with open(dst, 'w') as file:
        for item in name_list:
            if name_list.index(item)==0:
                file.write('{:5}\t'.format(item))
            elif name_list.index(item)==len(name_list)-1:
                file.write('{:15}\n'.format(item))
            else:
                file.write('{:15}\t'.format(item))
        for j in range(len(res_list)):
            item_list = res_list[j]
            for i in range(len(item_list)):
                item = item_list[i]
                if i==0:
                    file.write('{:5d}\t'.format(item))
                elif i==len(item_list)-1:
                    file.write('{:.9e}\n'.format(item))
                else:
                    file.write('{:.9e}\t'.format(item))
            
    print(res_list)
    

# if use_local_storage:
#     for runname in [runname_field,runname_lab]:
#         src = outdir + runname 
#         dst = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/' + runname 
        
#         if not os.path.exists(dst): 
#             shutil.copytree(src, dst)
#         else:
#             shutil.rmtree(dst)
#             shutil.copytree(src, dst)


# ... run postprocessing checks
shf.run_complete_check(runname_field, 
                      runname_lab, 
                      outdir, 
                      target_duration=ttot_field, 
                      include_duration_check=True, 
                      omit_saveSuff=True, 
                      omit_ipynb=True,
                     )

# ... compute profile data
cflx.prof_postproc_save(outdir, runname_field, runname_lab, postproc_prof_list)

# ... move to aws if this option is turned on
# [nothing happens if aws_save != 'move' or 'copy']
shf.to_aws(aws_save, 
           aws_bucket, 
           outdir, 
           runname_lab, 
           runname_field)
            