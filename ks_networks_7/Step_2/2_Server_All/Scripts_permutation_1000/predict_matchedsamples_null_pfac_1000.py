import glob
import h5py
import numpy as np
import os
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
import sys
import time
from scipy.stats import pearsonr

"""
This script will run ridge regression. To run it, you need:

output_dir: where you want your results saved, does not have to exist, we will make it
feature_path: path to the actual features, should be a subject x features npy array
phenotype_path: path to a csv with the to-be-predicted values
phenotype_name: str of name of column in phenotype_path you want
control_path: path to csv, regress these features from the features within each fold
fold_group: in this instance, this is the ABCD train and test split, A for one group, B for the other

"""

folds = 100
output_dir = sys.argv[1] #'homedir/results_matchedsamples_null_1000_mmddyyX/all_network/'
os.makedirs(output_dir,exist_ok=True)
feature_path = sys.argv[2]  #'homedir/tempfiles_matchedsamples_031324/features_all.npy'
phenotype_path = sys.argv[3] #'homedir/tempfiles_matchedsamples_031324/phenotypes_target.csv'
phenotype_name = sys.argv[4] #'General_PB1'
control_path = sys.argv[5] #'homedir/tempfiles_matchedsamples_031324/phenotypes_control.csv'
fold_group = 'matched_group' #sys.argv[6] # in this instance, this is the ABCD train and test split, A for one group, B for the other

"""
load the subject measures you want to control for
"""
phenotypes_control = pd.read_csv(control_path)
"""
load in all the feature weights
"""
features = np.load(feature_path).astype(np.float16)
"""
this is adapted from pennlinckit.utils.predict
"""
targets = pd.read_csv(phenotype_path)[phenotype_name].values.astype(np.float16)
fold_group = pd.read_csv(phenotype_path)[fold_group].values.astype(np.float16)
np.save('{0}/{1}_targets_null.npy'.format(output_dir,phenotype_name),targets)
assert targets.shape[0] == features.shape[0]

accuracy_testA = []
accuracy_testB = []
accuracy_boot =[]
prediction = np.zeros((targets.shape))

time = time.time() # set seed based on time to ensure random shuffling across 10x concurrent jobs
time = str(time)
time = time.split('.')
seed = int(time[1])
np.random.seed(seed)
print(seed)

for fold in np.arange(folds): #100 permutations
	# shuffle the outcome variable
	np.random.shuffle(targets)
	
	# Split up the data into train/test based on the fold group
	A = np.argwhere(fold_group==1) # <- save an index of everywhere the fold group is 1
	B = np.argwhere(fold_group==2) # <- save an index of everywhere the fold group is 2
	x_A = features[A]
	y_A= targets[A]
	x_B = features[B]
	y_B = targets[B]
	nuisance_A = phenotypes_control.values[A]
	nuisance_B = phenotypes_control.values[B]

	# Double check that the groups have no overlap
	assert np.intersect1d(A,B).size==0
	# Remove the unnecessary null dimension
	x_A = np.squeeze(x_A)
	y_A = np.squeeze(y_A)
	x_B = np.squeeze(x_B)
	y_B = np.squeeze(y_B)
	nuisance_A=np.squeeze(nuisance_A)
	nuisance_B=np.squeeze(nuisance_B)

	# First train on group A and test on group B:
	nuisance_model = LinearRegression() #make the nuisance model object
	nuisance_model.fit(nuisance_A,x_A) #fit the nuisance_model to training data
	x_A_resid = x_A - nuisance_model.predict(nuisance_A) #remove nuisance from training data
	x_B_resid = x_B - nuisance_model.predict(nuisance_B) #remove nuisance from test data
	m = RidgeCV(alphas=(1,10,100,500,1000,5000,10000,15000,20000)) #make the actual ridge model object, adding some super high reg strengths because we have so many features
	m.fit(x_A_resid,y_A) # fit the ridge model
	predicted_y_B = m.predict(x_B_resid) #apply the trained model to the test data
	coefs_B = m.coef_.astype(np.float16)

	np.save('{0}/{1}_prediction_testB_{2}.npy'.format(output_dir,phenotype_name,fold),predicted_y_B) # this is what will be used to calculate Haufe
	np.save('{0}/{1}_coefs_testB_{2}.npy'.format(output_dir,phenotype_name,fold),coefs_B)

	# Then train on group B and test on group A:
	nuisance_model = LinearRegression() #make the nuisance model object
	nuisance_model.fit(nuisance_B,x_B) #fit the nuisance_model to training data
	x_B_resid = x_B - nuisance_model.predict(nuisance_B) #remove nuisance from training data
	x_A_resid = x_A - nuisance_model.predict(nuisance_A) #remove nuisance from test data
	m = RidgeCV(alphas=(1,10,100,500,1000,5000,10000,15000,20000)) #make the actual ridge model object, adding some super high reg strengths because we have so many features
	m.fit(x_B_resid,y_B) # fit the ridge model
	predicted_y_A = m.predict(x_A_resid) #apply the trained model to the test data
	coefs_A = m.coef_.astype(np.float16)

	np.save('{0}/{1}_prediction_testA_{2}.npy'.format(output_dir,phenotype_name,fold),predicted_y_A) # this is what will be used to calculate Haufe
	np.save('{0}/{1}_coefs_testA_{2}.npy'.format(output_dir,phenotype_name,fold),coefs_A)

	# calculate prediction accuracy for A and B separately
	acc_testA = pearsonr(y_A, predicted_y_A)[0]
	acc_testB = pearsonr(y_B, predicted_y_B)[0]
	accuracy_testA.append(acc_testA)
	accuracy_testB.append(acc_testB)

	# Calculate the prediction accuracy for A and B combined
	acc_boot = pearsonr(np.concatenate((y_A,y_B)), np.concatenate((predicted_y_A,predicted_y_B)))[0] 
	accuracy_boot.append(acc_boot)

np.save('{0}/{1}_acc_null_testA.npy'.format(output_dir,phenotype_name),accuracy_testA) #save out prediction accuracy for A
np.save('{0}/{1}_acc_null_testB.npy'.format(output_dir,phenotype_name),accuracy_testB) #save out prediction accuracy for B
np.save('{0}/{1}_acc_boot_combined.npy'.format(output_dir,phenotype_name),accuracy_boot) #save out prediction accuracy for A and B combined



