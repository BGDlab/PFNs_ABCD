import glob
import h5py
import numpy as np
import os
import pandas as pd
from sklearn import preprocessing
import sys

"""
hi. this is going to read these in, you basically just need 
subjects: path to subjects!
data_path: path to the actual features
phenotypes: path to a csv with the to-be-predicted and the noise

"""
data_path ='/cbica/projects/PFN_ABCD/input_files/{0}/IndividualParcel_Final_sbj1_comp17_alphaS21_1_alphaL300_vxInfo1_ard0_eta0/final_UV.mat' # feature path
control =['sex','meanFD','abcd_site','interview_age'] # covariate variable names
try: network = int(sys.argv[1])
except: network = sys.argv[1]
outdir = sys.argv[2]   #tmpdir
os.makedirs(outdir,exist_ok=True)
"""
load the subject measures 
"""
phenotypes = pd.read_csv('./data_for_ridge_PRS.csv') #load covariate and outcome data
#phenotypes['subjectkey'] = np.core.defchararray.strip(phenotypes.subjectkey.values.astype(str),'NDAR_') #strip
phenotypes = phenotypes.sort_values('subjectkey') #sort 
phenotypes['meanFD'] = phenotypes['meanFD'].values.astype(np.float16) #saves space

"""
convert categorical vars to numerical, ie dummy encoding, for sites and sex
"""
for c in ['sex','abcd_site']:
    dumb = pd.get_dummies(phenotypes[c]).astype(np.int8)
    phenotypes = pd.merge(left=phenotypes,right=dumb,left_index=True,right_index=True)
    for d in dumb.columns:
        control.append(d) #we get a lot more columns so we put them in
    control.remove(c) # remove the original categorical

"""
load in all the feature weights
"""
features = []
for subject in phenotypes.subjectkey.values:
    if os.path.exists(data_path.format(subject)) == False:
        phenotypes = phenotypes.drop(phenotypes.index[phenotypes.subjectkey==subject])
        continue
    f = h5py.File(data_path.format(subject), 'r')
    f = np.array(f[f['V'][0][0]]).astype(np.float16)
    if type(network) == int:
        f = f[network]
    features.append(f)
features = np.array(features)
features = features.reshape(features.shape[0],-1) #subject by flat array shape to be flat for sklearn

np.save('/{0}/features_{1}.npy'.format(outdir,network),features.astype(np.float16))



phenotypes[control].to_csv('{0}/phenotypes_control.csv'.format(outdir))
phenotypes[['PRS_1','PRS_2','rel_family_id','matched_group']].to_csv('{0}/phenotypes_target.csv'.format(outdir))
