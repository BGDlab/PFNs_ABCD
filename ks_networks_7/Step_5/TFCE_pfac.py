import os
import numpy as np
import nibabel as nb
from scipy.sparse import lil_matrix
import string
import joblib
from tqdm import tqdm

#LOAD IN ACTUAL WEIGHT MAP AVG AND SUM
# observed weight map
weights_act_A = np.load("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_5/BASELINE/Weight_Haufe_031324/General_PB1_A_Weight_haufetrans_all.npy")
weights_act_B = np.load("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_5/BASELINE/Weight_Haufe_031324/General_PB1_B_Weight_haufetrans_all.npy") 
weights_act_avg = (weights_act_A + weights_act_B) / 2
weights_act_avg_abs = np.abs(weights_act_avg) #abs value
weights_act_avg_abs = weights_act_avg_abs.reshape(17,59412).T #reshape to network structure
weights_act_sum = weights_act_avg_abs.sum(axis=1) #sum across network
np.savetxt("weight_map_HaufeAAB/sig_testing/pfac_vertex_SOW_python.csv", weights_act_sum, delimiter=",", header="SOW", comments="")
print(f"Loaded actual summed weights")
print("max data:", weights_act_sum.max())

#LOAD NULL WEIGHT AVG AND SUM
#Path to null weight directories
base_dir = "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_5/BASELINE/Haufe_Nulls_031824/"
# Pre‐allocate list to collect each permutation’s weight vector
null_weights_list = []
# Loop through directories A–J
for dir_letter in string.ascii_uppercase[:10]:  # 'A' through 'J'
    dir_path = os.path.join(base_dir, f"Weight_Haufe_Nulls_031824_{dir_letter}")
    if not os.path.isdir(dir_path):
        raise FileNotFoundError(f"Directory not found: {dir_path}")
    
    # Loop through files 0.npy … 99.npy
    for i in range(0, 100):
        fn_A = f"General_PB1_A_Weight_haufetrans_all_{i}.npy"
        fp_A = os.path.join(dir_path, fn_A)
        fn_B = f"General_PB1_B_Weight_haufetrans_all_{i}.npy"
        fp_B = os.path.join(dir_path, fn_B)
        if not os.path.isfile(fp_A):
            raise FileNotFoundError(f"File A not found: {fp_A}")
        if not os.path.isfile(fp_B):
            raise FileNotFoundError(f"File B not found: {fp_B}")
            
        # Load and append
        null_weights_A = np.load(fp_A, mmap_mode="r")          # shape (59412,)
        null_weights_B = np.load(fp_B, mmap_mode="r")          # shape (59412,)
        null_weights_avg = (null_weights_A + null_weights_B) / 2
        null_weights_avg_abs = np.abs(null_weights_avg) #abs value
        null_weights_avg_abs = null_weights_avg_abs.reshape(17,59412).T #reshape to network structure
        null_weights_sum = null_weights_avg_abs.sum(axis=1) #sum across networks
        null_weights_list.append(null_weights_sum)
# Stack into a single array of shape (1000, 59412)
null_weights = np.stack(null_weights_list, axis=0)
print(f"Loaded null weights: {null_weights.shape[0]} permutations, each of length {null_weights.shape[1]}")

#LOAD NON-MEDIAL WALL MASK
n_mask = weights_act_sum.shape[0]
template = nb.load("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Step_6/hardparcel_group.dscalar.nii")
left_inds = []
right_inds = []
for model in template.header.get_index_map(1).brain_models:
    arr = np.asarray(model.vertex_indices)
    if model.brain_structure == "CIFTI_STRUCTURE_CORTEX_LEFT":
        left_inds = arr
    elif model.brain_structure == "CIFTI_STRUCTURE_CORTEX_RIGHT":
        right_inds = arr
nverts_h = 32492
right_inds = right_inds + nverts_h
mask_inds = np.concatenate((left_inds, right_inds))
assert mask_inds.shape[0] == n_mask #should be 59,412
assert weights_act_sum.shape[0] == mask_inds.shape[0]
assert np.unique(mask_inds).size == mask_inds.size #checking for duplicate indices

#LOAD FSLR SURFACE TOPOLOGY
# Load left hemisphere
surf_l = nb.load("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii")
verts_l = surf_l.darrays[0].data      # shape (32492, 3)
faces_l = surf_l.darrays[1].data.astype(int)  # shape (64980, 3)
# Load right hemisphere
surf_r = nb.load("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Satterthwaite/ks_networks_7/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
verts_r = surf_r.darrays[0].data      # shape (32492, 3)
faces_r = surf_r.darrays[1].data.astype(int)  # shape (64980, 3)
# Combine into a single mesh
# Offset right-face indices by the number of left vertices
nverts_l = verts_l.shape[0] #32492
faces_r_global = faces_r + nverts_l
# Stack vertices and faces
verts = np.vstack((verts_l, verts_r))
faces = np.vstack((faces_l, faces_r_global))

#BUILD ADJACENCY BASED ON SURFACE TOPOLOGY
def build_connectivity(faces, nverts):
    adj = lil_matrix((nverts, nverts), dtype=bool)
    for i, j, k in faces:
        adj[i, j] = adj[j, i] = True
        adj[i, k] = adj[k, i] = True
        adj[j, k] = adj[k, j] = True
    return adj.tocsr()

adj_full = build_connectivity(faces, nverts=verts.shape[0])
#restrict adjacency to non-medial wall mask
adj = adj_full[mask_inds, :][:, mask_inds]
assert adj.shape == (n_mask, n_mask)

#TFCE function
def tfce(data, connectivity, E=0.5, H=2.0, dh=0.01):
    """
    Pure-Python TFCE without nimare.
    `data`: 1D array of length n_mask
    `connectivity`: sparse adjacency (csr_matrix of shape n_mask x n_mask)
    """
    n = data.size #n of vertices
    tfce_score = np.zeros(n, float) #TFCE score starts at zero
    # thresholds from dh up to max(data)
    thresholds = np.arange(dh, data.max() + dh, dh) #building list of thresholds, dh is the increment of steps up to the max data value
    # pre-extract neighbor lists
    neighbors = [connectivity[i].indices for i in range(n)]
    for h in thresholds: #loop over all thresholds
        mask = data >= h #mask of vertices that exceed the current threshold
        visited = np.zeros(n, bool) #track already visited vertices
        for i0 in np.where(mask & ~visited)[0]: #find unvisited vertices, then grow clusters
            stack, cluster = [i0], []
            while stack:
                v = stack.pop()
                if not visited[v] and mask[v]:
                    visited[v] = True
                    cluster.append(v)
                    stack.extend(neighbors[v])
            if cluster:
                size = len(cluster)
                tfce_score[cluster] += (h**H) * (size**E) * dh
    return tfce_score

#Actual TFCE computation
# H and E are the TFCE exponents (defaults H=2.0, E=0.5)
# dh is the step size in your threshold integral, set to 0.01
tfce_act = tfce(
    weights_act_sum,
    connectivity=adj
)  # shape (n_vertices,)
assert tfce_act.shape[0] == n_mask
np.savetxt("weight_map_HaufeAAB/sig_testing/pfac_TFCE_act.csv", tfce_act, delimiter=",", header="tfce_score", comments="")

#Null TFCE computation loop w/ parallelization
n_perm = null_weights.shape[0] #1000
def max_tfce(w):
    t = tfce(w, connectivity=adj)
    return np.max(np.abs(t))

max_tfce_null = joblib.Parallel(n_jobs=8)(
    joblib.delayed(max_tfce)(null_weights[i])
    for i in tqdm(range(n_perm), desc="TFCE nulls")
)
max_tfce_null = np.array(max_tfce_null)


# include a “+1” in numerator & denominator for a proper permutation p-value
p_fwe = (1 + np.sum(np.abs(tfce_act)[None, :] <= max_tfce_null[:, None], axis=0)) / (1 + n_perm)

# p_fwe is your 1-D array of length n_mask
np.savetxt("weight_map_HaufeAAB/sig_testing/pfac_TFCE_p_vals.csv", p_fwe, delimiter=",", header="p_fwe", comments="")
