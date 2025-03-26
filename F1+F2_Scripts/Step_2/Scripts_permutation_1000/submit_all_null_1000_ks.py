import os
import sys
import time
import numpy as np
"""
syntax:  python submit_all_null_1000_ks.py <mmddyyX> X=letter from A-J (10 jobs)

*****  this is only run network all null  *****

this is the only script you actually should call, it will submit a job
this job makes the features, and then does the regression. 
it assigns the current running dir to home_dir variable. It refers to/creates 
   tmp dir = <current dir>/tempfiles_matchedsamples_040524
   result dir = <current dir>/results_matchedsamples_null_1000_<mmddyyX> 
"""

dirdate = sys.argv[1]   # the command line parameter as mmddyyX 
homedir = os.getcwd()   #get the current running dir '/cbica/projects/PFN_ABCD/scripts/ridge/Step_2/ridge_matchedsamples' #sys.argv[1] 

tmpdir = '{0}/tempfiles_matchedsamples_040524'.format(homedir)
sge_dir = '{0}/sge_{1}'.format(homedir,dirdate)  # dir for the processors std output and system error output files

os.makedirs(sge_dir,exist_ok=True)  # Create the output dir if it doesn't exist

networks = np.array(['all']) 
network = 'all'

GB = '400G'

for outcome in ['PRS_1','PRS_2']:
    print('qsub -l h_vmem={0},s_vmem={0} -N {1}{4} -R y -V -j y -b y -o ./sge_{6}/ -e ./sge_{6}/ python {2}/ks_proc_predict_null_1000.py {2} {3} {4} {5} {6}'.format(GB,outcome[1],homedir,tmpdir,network,outcome,dirdate))
    os.system('qsub -l h_vmem={0},s_vmem={0} -N {1}{4} -R y -V -j y -b y -o ./sge_{6}/ -e ./sge_{6}/ python {2}/ks_proc_predict_null_1000.py {2} {3} {4} {5} {6}'.format(GB,outcome[1],homedir,tmpdir,network,outcome,dirdate))
    #time.sleep(3)
