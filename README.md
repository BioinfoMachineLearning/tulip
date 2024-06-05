## TULIP
The template-based modeling for accurate ligand-protein complex structure prediction (TULIP).

``Input`` : Query protein, Ligand molecule ( .SMILES ), Protein Template Hits

``Output`` : ``.sdf`` file of ligand molecule that potentially binds to input query protein

## Setup Environment
1. The required packages are found in ``environment.yml`` 

``conda env create -f environment.yml``

2. Once the installation is completed, activate the conda environment

``conda activate tulip``

3. The pipeline makes use of UCSF Chimera and PyRosetta. Please install them through their websites

[UCSF Chimera Download](https://www.cgl.ucsf.edu/chimera/download.html) , this link will take users to: ```https://www.cgl.ucsf.edu/chimera/download.html```

[PyRosetta Download](https://www.pyrosetta.org/downloads) , this link will take users to: ```https://www.pyrosetta.org/downloads```


To adjust the target ligand's binding pose and orientation by rotation and translation, TULIP uses LS-align to align the target ligand with template ligands of higher similarity by both flexible and rigid alignments.


The LS-align tool can be downloaded and compiled into local machine from the below link:

[LS-align](https://zhanggroup.org/LS-align/) , this link will take users to: ```https://zhanggroup.org/LS-align/```

## Run

To run TULIP pipeline run the following:

``python3 complex_generator.py --config=../configs/configs.yml``


``configs/configs.yml`` contains the configuration of the pipeline. Enter the correct paths in config file.

## Post Processing

Once the ligand templates are identified using the above script, use the script below to run the post-processing for multi-atom ligands. Update the ``ls_align_path`` variable in the script to the path where LS-align is installed.

``python3 post_process.py ``

For ions and small molecules run:

``python3 process_small_mol.py ``

## Clustering 
To perform clustering in case we have same ligand SMILES for a protein.

``python3 cluster_files.py``

## Outputs
TULIP returns the identified ligand molecules in ranks. Example : `rank_0.sdf`, `rank_1.sdf`, and so on for each protein.



