
## TULIP
Template-based modeling for accurate ligand-protein complex structure prediction (TULIP) is developed by BML lab at University of Missouri, Columbia for protein-ligand complex modeling in CASP15.

## Setup Environment
1. The required packages are found in ``environment.yml`` 

``conda env create -f environment.yml``

2. Once the installation is completed, activate the conda environment

``conda activate tulip``

3. The pipeline makes use of UCSF Chimera and PyRosetta. Please install them through their websites

[UCSF Chimera Download](https://www.cgl.ucsf.edu/chimera/download.html) , this link will take users to: ```https://www.cgl.ucsf.edu/chimera/download.html```

[PyRosetta Download](https://www.pyrosetta.org/downloads) , this link will take users to: ```https://www.pyrosetta.org/downloads```

## Run
To run the whole pipeline run the following:

    python3 complex_generator.py --config=../configs/run.yml 

``configs/run.yml`` contains the configuration of the pipeline, example: change ``input_dir`` and ``output_dir`` from ``run.yml`` file.

The outputs are currently saved in 
``'./outputs'``

## Input directory
The data are organized in the following way in the input directory:

```

data/inputs
└───Target 1                        # target protein-ligand complex name
    |   predicted                   
        └───query.pdb               # predicted receptor protein structure using MULTICOM or MULTICOM_qa
    |   templates
        └───pdb1ctn.ent             # identified templates using Foldseek
        └───pdb1d2k.ent
        └───...
        └───...
        └─── pdb7vrg.ent
    |   evalue.m8                   # Estimated E-values for each match from Foldseek
    |   lig.smiles                  # Ligand SMILES, from CASP15 website
└───Target 2
    |   predicted                   
        └───query.pdb               
    |   templates
        └───pdb1ctn.ent             
        └───pdb1d2k.ent
        └───...
        └───...
        └─── pdb7vrg.ent
    |   evalue.m8                   
    |   lig.smiles 
    
...
```


## Output directory
``'data/outputs/Target1'`` 
there are files which are generated during runtime of the program. Some main file name convention are as
follows:
1. ``ligand_template_2_0.pdb`` -> ligand extracted from PDB template `2` and is in coordinate space of predicted structure `0`
2. ``PEE_2_0.pdb`` -> Single pdb file for unique ligand `PEE` which is identified from ``ligand_template_2_0.pdb``. This is important because ligands extracted from PDB template can have multiple ligands in them, so we separate them and save them in different file with ligand name as prefix. here, in the filename, the number `2` refers to the template `2` and `0` is predicted structure `0`.
3. ``PEE_mol_sim_Finger_0.82_Dice_0.829_2_0_smiles_0.mol`` -> final ligand `.mol` file in coordinate space of predicted structure `0` using template `2` with Tanimoto similarity of `0.82` and Dice similarity of `0.829` for smiles string `0`.

## LS-align:
To adjust the target ligand's binding pose and orientation by rotation and translation, TULIP uses LS-align to align the target ligand with template ligands of higher similarity by both flexible and rigid alignments.

**Downloading the tool: Locally**

The LS-align tool can be downloaded and compiled into local machine from the below link:

[LS-align](https://zhanggroup.org/LS-align/) , this link will take users to: ```https://zhanggroup.org/LS-align/```

<ins>**Running LSAlign tool: Locally compiled**</ins>

Once, the tool is complied, the tool can be run as follows:

``src/LSalign QUERYs.mol2 TEMPLs.mol2 -rf 1 -o flexible.mol2 `` (runs the flexible alignment)

`` src/LSalign QUERYs.mol2 TEMPLs.mol2 -rf 0 -o rigid.mol2`` (runs the rigid alignment)

Here, ``QUERYs.mol2`` is the target ligand molecule file and ``TEMPLs`` is the template ligand molecule file. The rigid and flexible flag is toggled by ``-rf``, where flag``1`` runs flexible alignment and flag ``0`` (default) runs rigid alignment. ``-o`` is the output location and name of the output.


<ins>**Running LSAlign tool : Web based**</ins>

To run the web based tool of LSAlign, please visit the website below and upload target and template ligand structures.

[LS-align](https://zhanggroup.org/LS-align/) , this link will take users to: ```https://zhanggroup.org/LS-align/```


## Fitted Ligand for Target T1158v3 using TULIP

![T1158v3](example/T1158v3_Human.gif)
