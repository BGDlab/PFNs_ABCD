import os
import sys
import time
import numpy as np
"""
syntex:  python submit_all_ks_scanner.py <mmddyy>

*****  this is only run network all (submit_ks.py runs network all and 0 to 16 (network 1-17))  *****

this is the only script you actually should call, it will submit a job
this job makes the features, and then does the regression. 
it assigns the current running dir to home_dir variable. It creates 
   tmp dir = <current dir>/tempfiles_matchedsamples_<mmddyy>
   result dir = <current dir>/results_matchedsamples_<mmddyy> 
"""

dirdate = sys.argv[1]   # the command line parameter as mmddyy 
homedir = os.getcwd()   #get the current running dir '/cbica/projects/PFN_ABCD/scripts/ridge/Step_2/ridge_matchedsamples' #sys.argv[1] 

tmpdir = '{0}/tempfiles_matchedsamples_{1}'.format(homedir,dirdate)  #It will be created if not exist    (or os.environ['TMPDIR'])
slurm_dir = '{0}/slurm_{1}'.format(homedir,dirdate)  # dir for the processors std output and error files

os.makedirs(slurm_dir, exist_ok=True)  # Create the output dir if it doesn't exist

networks = np.array(['all']) 
network = 'all'

# Run the preprocessing step
print('python {0}/preprocess_ks_scanner.py {1} {2}'.format(homedir,network,tmpdir))
os.system('python {0}/preprocess_ks_scanner.py {1} {2}'.format(homedir,network,tmpdir))

GB = '200G'

for outcome in ['General_PB1']:
    job_name = f'{outcome[:11]}_{network}'
    slurm_command = (
        f'sbatch --mem={GB} --job-name={job_name} --output={slurm_dir}/%x_%j.out --error={slurm_dir}/%x_%j.err '
        f'--wrap="python {homedir}/ks_proc_predict_scanner.py {homedir} {tmpdir} {network} {outcome} {dirdate}"'
    )
    
    print(slurm_command)
    os.system(slurm_command)
    # Optionally, add a small delay between job submissions
    # time.sleep(3)
