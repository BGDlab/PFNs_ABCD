import sys
import os

homedir = sys.argv[1]  #'/cbica/projects/PFN_ABCD/scripts/ridge/Step_2/ridge_matchedsamples'
tmpdir = sys.argv[2]   #'/cbica/projects/PFN_ABCD/scripts/ridge/Step_2/ridge_matchedsamples/tempfiles_matchedsamples_mmddyyX'
network = sys.argv[3]  # all
phenotype_name = sys.argv[4]   # PRS
dirdate = sys.argv[5]   # the date as mmddyyX

output_dir = homedir + '/results_matchedsamples_null_1000_{0}/{1}_network/'.format(dirdate,network)
feature_path = '{0}/features_{1}.npy'.format(tmpdir,network)
phenotype_path = '{0}/phenotypes_target.csv'.format(tmpdir)
control_path = '{0}/phenotypes_control.csv'.format(tmpdir)

"""
this script does the actual ridge regression
"""
print('python {0}/predict_matchedsamples_null_1000.py {1} {2} {3} {4} {5} rel_family_id'.format(homedir,output_dir,feature_path,phenotype_path,phenotype_name,control_path))
os.system('python {0}/predict_matchedsamples_null_1000.py {1} {2} {3} {4} {5} rel_family_id'.format(homedir,output_dir,feature_path,phenotype_path,phenotype_name,control_path))
