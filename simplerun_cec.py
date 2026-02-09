import json
import os,shutil,sys
from pathlib import Path

import fsspec
import numpy as np
import spinup,get_int_prof,get_soilpH_time

# --- read in helper functions from ew-workflows
from ew_workflows import scepter_helperFxns as shf
from ew_workflows import cflx_proc as cflx

def read_json(path):
    """
    Read a JSON file from local disk or S3.
    
    Args:
        path (str or Path): Local path or S3 path (s3://bucket/key)
    
    Returns:
        dict: The loaded JSON data.
    """
    path = str(path)  # ensure string
    if path.startswith('s3://'):
        fs = fsspec.filesystem('s3')
        with fs.open(path, 'r') as f:
            data = json.load(f)
    else:
        with open(path, 'r') as f:
            data = json.load(f)
    return data

        
def simplerun(**rundict):

    # --- populate inputs
    # 
    # outdir_src = rundict.get('outdir', '/home/jovyan/SCEPTER/scepter_output')
    # use_local_storage = False # (True for GT cluster)
    rstrt = rundict.get('rstrt', 'self')

    # (keeping from Yoshi's original script)
    # ---- from Tipping_Hurley.dat ---- 
    # Na+ + X- = NaX
    # log_k           0.0
    #
    # K+ + X- = KX
    # log_k           0.7
    #
    # H+ + X- = HX
    # log_k           1.0
    #
    # Ca+2 + 2X- = CaX2
    # log_k           0.8
    #
    # Mg+2 + 2X- = MgX2
    # log_k           0.6
    #
    # Al+3 + 3X- = AlX3
    # log_k           0.67
    # ---------------------------------
    
    # ---- how SCEPTER parameterizes ---- 
    # X-Na + H+ = Na+(aq) + X-H (logkhna)
    # where 
    #       logkhna   = logkhna_0 * gamma
    #       gamma = 10^(  alpha * f[X-H] ) 
    # 
    # assume database is defined with alpha = 0
    # logkhna   = 1.0 - 0.0 = 1.0
    # logkhk    = 1.0 - 0.7 = 0.3
    # 
    # multi valent cations 
    # X2-Ca + 2H+ = Ca++(aq) + 2X-H 
    # 
    # logkhca   = (1.0 - 0.8)*2 = 0.4
    # logkhmg   = (1.0 - 0.6)*2 = 0.8
    # logkhal   = (1.0 - 0.67)*3 = 0.99
    # ---------------------------------- 
    
    # (keeping Yoshi's notes below on pore volume exchange math)
    # 0.0011 mol/1kgw of exchanger 
    # assume poro m3/m3 porosity, sat as m3/m3 water saturation, 1 - poro m3/m3 of solid phase
    # poro * sat * 1,000 kg/m3 of solution assuming 1,000 kg/m3 = 1 g/cm3 solution density
    # (1-poro) * rho kg/m3 of solid given the density rho in kg/solid kg
    # assuming cec in units of cmol/kg, total exchange sites in soil is (1-poro) * rho * cec cmol/m3
    # exchange sites for a given solution mass is then given as 
    #   (1-poro) * rho * cec * 0.01 / ( poro * sat * 1,000 )  mol/kgw 
    # and this must be equal to 0.0011 mol/1kgw
    #   (1-poro) * rho * cec * 0.01 / ( poro * sat * 1,000 ) = 0.0011
    # [ yoshi's cec calc depends on rho to align with PHREEQC benchmark ]
    # rho = 258.162/99.52 * 1000 # kg/m3
    # cec = 0.0011 * (0.5 * 1 * 1000 )/ ( (1.-0.5) * rho * 0.01 )
    
    
    if rstrt != 'self':
        print("TKTK write code to check if restart is on s3 and move it to local !!")
        return "Need to set up code to handle restart from spinup for `simplerun_cec.py` (just move it from s3 to local)"
    
    
    spinup.run_a_scepter_run(**rundict)
   
def main():
    # check that a path was provided
    if len(sys.argv) < 2:
        print("Usage: python3 simplerun.py path/to/rundict.json")
        sys.exit(1)

    # read in the run dict
    # --- read in the rundict (assume json)
    rundict_path = sys.argv[1]  # get the first command-line argument
    rundict = read_json(rundict_path)
    
    # change dir for relative path integration (esp. in make_inputs..)
    os.chdir(rundict['model-dir'])
    
    # run
    simplerun(**rundict)

    # postproc
    outdir = rundict.get('outdir_src', "NOOUTDIR!!")
    runname = rundict.get('runname', 'NORUNNAME!!')
    # ... compute profile data
    postproc_prof_list = rundict.get('postproc_prof_list', 'all')
    cflx.prof_postproc_save(outdir, runname, 
                            runname_lab=None, postproc_prof_list=postproc_prof_list)
    
    # ... move to aws if this option is turned on
    # [nothing happens if aws_save != 'move' or 'copy']
    aws_save = rundict.get('aws_save', None)
    aws_bucket = rundict.get('aws_bucket', None)
    shf.to_aws(aws_save, 
               aws_bucket, 
               outdir, 
               runname_lab=None, 
               runname_field=runname)

   
if __name__ == '__main__':
    main()
    
    