"""
Runs clustering in case we have same ligand SMILES for a protein.

"""


import os
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.cluster import AgglomerativeClustering
import numpy as np
from scipy.spatial.distance import cdist
import shutil


num_clusters = 2
protein = "T1187" # protein name
ligand_name = "001" # ligand id
file_dir = f"/outputs/{protein}/{ligand_name}" 
des = f"/outputs/{protein}/final"
sdf_files = [os.path.join(file_dir, file) for file in os.listdir(file_dir) if file.endswith('.sdf')]

# Step 1: Read .sdf files using RDKit
ligands = []
for filename in sdf_files:
    suppl = Chem.SDMolSupplier(filename)
    for mol in suppl:
        if mol:
            ligands.append(mol)

# Step 3: Compute pairwise distances between ligands
num_ligands = len(ligands)
coords = np.array([lig.GetConformer().GetPositions() for lig in ligands])
distances = cdist(coords.reshape(num_ligands, -1), coords.reshape(num_ligands, -1))

# Step 4: Generate clusters
clustering = AgglomerativeClustering(n_clusters=num_clusters, affinity='precomputed', linkage='complete').fit(distances)
labels = clustering.labels_

# Step 5: Identify ligands in clusters
clusters = {}
for i, label in enumerate(labels):
    if label not in clusters:
        clusters[label] = []
    clusters[label].append(sdf_files[i])

# Print ligands in each cluster
for cluster_id, ligands_in_cluster in clusters.items():
    print(f'Cluster {cluster_id + 1}:')
    for ligand_file in ligands_in_cluster:
        print(os.path.basename(ligand_file))
        
cluster_dict = {}

for cluster_id, ligands_in_cluster in clusters.items():
    sorted_ligands = sorted([os.path.basename(ligand_file) for ligand_file in ligands_in_cluster], key=lambda x: int(x.split('_')[1].split('.')[0]))
    cluster_dict[cluster_id + 1] = sorted_ligands

print(cluster_dict)

for k,v in cluster_dict.items():
    for li in range(len(v)):
        src_file = f"{file_dir}/{v[li]}"
        des_file = f"{des}/{ligand_name}_{k}"
        des_file1 = f"{des}/{ligand_name}_{k}/rank_{li}.sdf"

        if not os.path.exists(des_file):
            os.makedirs(des_file)
        shutil.copy(src_file, des_file1)
