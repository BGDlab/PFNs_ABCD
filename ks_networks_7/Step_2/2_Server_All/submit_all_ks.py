import os
import sys
import time
import numpy as np
"""
syntex:  python submit_all_ks.py <mmddyy>

*****  this is only run network all (submit_ks.py runs network all and 0 to 16 (network 1-17))  *****

this is the only script you actually should call, it will submit a job
this job makes the features, and then does the regression. 
it assigns the current running dir to home_dir variable. It creates 
   tmp dir = <current dir>/tempfiles_matchedsamples_<mmddyy>
   result dir = <current dir>/results_matchedsamples_<mmddyy> 
the network final_UV.mat files path need to be updated in preprocess_ks.py code
thing else should run! 
"""

dirdate = sys.argv[1]   # the command line parameter as mmddyy 
homedir = os.getcwd()   #get the current running dir '/cbica/projects/PFN_ABCD/scripts/ridge/Step_2/ridge_matchedsamples' #sys.argv[1] 

tmpdir = '{0}/tempfiles_matchedsamples_{1}'.format(homedir,dirdate)  #It will be created if not exist    (or os.environ['TMPDIR'])
sge_dir = '{0}/sge_{1}'.format(homedir,dirdate)  # dir for the processors std output and system error output files

os.makedirs(sge_dir,exist_ok=True)  # Create the output dir if it doesn't exist

networks = np.array(['all']) 
network = 'all'


# following call preprocess_ks.py is moved from ks_proc_predict.py to avoid repeating for each psy_score run which causes the processor to die due to system resource issues.
print('python {0}/preprocess_ks.py {1} {2}'.format(homedir,network,tmpdir))
os.system('python {0}/preprocess_ks.py {1} {2}'.format(homedir,network,tmpdir))

GB = '200G'

for psy_score in ['General_PB1']:
    print('qsub -l h_vmem={0},s_vmem={0} -N {1}{4} -R y -V -j y -b y -o ./sge_{6}/ -e ./sge_{6}/ python {2}/ks_proc_predict.py {2} {3} {4} {5} {6}'.format(GB,psy_score.split('_')[1],homedir,tmpdir,network,psy_score,dirdate))
    os.system('qsub -l h_vmem={0},s_vmem={0} -N {1}{4} -R y -V -j y -b y -o ./sge_{6}/ -e ./sge_{6}/ python {2}/ks_proc_predict.py {2} {3} {4} {5} {6}'.format(GB,psy_score.split('_')[1],homedir,tmpdir,network,psy_score,dirdate))
    #time.sleep(3)
