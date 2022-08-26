"""
author: nabin 
timestamp: Wed May 04 2022 2:02 AM

- makes mol file from smiles -> for equibind

"""
import subprocess

import os
import os.path as osp
from rdkit import Chem
from rdkit.Chem import AllChem


def read_smiles(casp_ligand):
    smiles_found = list()
    ligand_file = open(casp_ligand)
    ligand_smiles = ligand_file.readlines()
    for line in ligand_smiles:
        line_split = line.split("   ")
        print(line_split)
        lig_smi = line_split[1]
        lig_smi = lig_smi.split(" ")
        lig_name = lig_smi[0]
        smi = lig_smi[1]
        print(f"CASP Ligand Name: {lig_name}")
        # if lig_name == "LIG":
        smiles_found.append(smi.strip("\n"))
    return smiles_found


def ligands_mol_from_smiles(casp_ligand):
    # smile = [ligand for ligand in ligand_path if ligand.endswith('.smiles')][0]
    smiles_list = read_smiles(casp_ligand)
    print(smiles_list)
    for smi in smiles_list:
        parent = os.path.join(casp_ligand, os.pardir)
        save_path = os.path.join(os.path.abspath(parent),
                                 f"equi_lig_{smiles_list.index(smi)}.mol")

        smi_mol = Chem.MolFromSmiles(smi)
        smi_1 = AllChem.EmbedMolecule(smi_mol)

        print(Chem.MolToMolBlock(smi_mol), file=open(save_path, 'w+'))

        # with open(osp.join(osp.abspath(parent), 'smiles_from_casp_for_equi.smiles'), 'w') as f:
            # f.write(smi)
        obabel_scripts = open(osp.join(osp.abspath(parent), 'obabel_script.cmd'), 'w')

        save_sdf = save_path.split(".")[0]
        obabel_scripts.write('obabel ' + '-imol ' + save_path + ' -osdf ' + f'-O{save_sdf}.sdf')
        obabel_scripts.close()

        command_complete = False
        while not command_complete:
            try:
                subprocess.run(['chmod', 'u+x', osp.join(osp.abspath(parent), 'obabel_script.cmd')])
                subprocess.run(osp.join(osp.abspath(parent), 'obabel_script.cmd'), shell=True)
                command_complete = True
            except FileNotFoundError as error:
                raise error
        os.remove(osp.join(osp.abspath(parent), 'obabel_script.cmd'))
        # os.remove(osp.join(osp.abspath(parent), 'smiles_from_casp_for_equi.smiles'))



