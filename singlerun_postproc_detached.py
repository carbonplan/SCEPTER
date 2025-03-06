# %% 
# -----------------------------------------
# 
# Run scepter in a scenario where the 
# run directory already exists
# (unlike the other scripts, this one does 
#  not create a new dir and populate it.. 
#  it just runs scepter and postprocessing
#  for a given dir)
# 
# -----------------------------------------
import os 
import re
import sys

sys.path.append(os.path.abspath('/home/tykukla/aglime-swap-cdr/scepter/setup'))
# import module
import scepter_helperFxns as shf
import cflx_proc as cflx


outdir = "/home/tykukla/SCEPTER/scepter_output/"
# runname = "_cectest_cec41"
# runname = "_cectest_alph20"
runname = "lowFert_cc_basev3_site_311a_app_30p0_psize_125_cc_field_tau15p0"


# %% 
# --- some inputs for postproc
fdust, fdust2 = 60, 6
dustsp, dustsp_2nd = 'cc', 'amnt'
postproc_prof_list = ["adsorbed_percCEC", "bulksoil", "aqueous", "gas", "solid_weightPercent"]
runname_lab = None

# %% 
# --- run scepter 
print('starting a scepter run...')
os.system("chmod +x " + os.path.join(outdir,runname,'scepter'))  # grant permissions
os.system(os.path.join(outdir,runname,'scepter'))


# %% 
# --- run postprocessing
runname_field = runname
# ... run postprocessing checks
shf.run_complete_check(runname_field, 
                      runname_lab, 
                      outdir, 
                      target_duration=15, 
                      include_duration_check=True, 
                      omit_saveSuff=True, 
                      omit_ipynb=True,
                     )

# ... update the dust flux file so it's usable without needing input from frame.in
shf.dustflx_calc(outdir, runname, fdust, fdust2, dustsp, dustsp_2nd)

# ... compute cdr-relevant fluxes
cflx.cflx_calc(outdir, runname, [dustsp, dustsp_2nd])

# ... compute profile data
cflx.prof_postproc_save(outdir, runname, runname_lab, postproc_prof_list)

# ... move to aws if this option is turned on
# [nothing happens if aws_save != 'move' or 'copy']
# shf.to_aws(aws_save, 
#            aws_bucket, 
#            outdir, 
#            runname_lab, 
#            runname_field)



# %%
