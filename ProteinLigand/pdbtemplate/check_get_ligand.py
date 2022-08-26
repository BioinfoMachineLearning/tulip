"""
author: nabin 
timestamp: Wed Apr 20 2022 3:50 PM

- check presence of ligand in the PDB template file and return True if ligand is present
"""
import os
import os.path as osp
from pdbtemplate import save_pdb_format as save_pdb

from pyrosetta import *

init()


def get_save_ligands(pdb_template, save_path):
    # os.makedirs(save_path, exist_ok=True)
    count = 0
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Extracting Ligands")
    pose_pdb = pose_from_pdb(pdb_template)
    print(
        f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>The total number of residue found for the pdb template {pdb_template} protein is : {pose_pdb.total_residue()}")
    # save_pdb.save_titles(file_path=save_path, target=1, remark="used PDB template",
    # method="identify ligands from template", model_number=1)
    for i in range(1, pose_pdb.total_residue() + 1):
        if pose_pdb.residue(i).is_ligand():
            print(f"The identified ligand outside is :{pose_pdb.residue(i).name()}")
        if pose_pdb.residue(i).is_ligand() and len(pose_pdb.residue(i).name()) <= 7:
            # if pose_pdb.residue(i).name()[-3:] == casp_ligand:
            print(f"The identified ligand inside is :{pose_pdb.residue(i).name()}")
            print(f"The identified chain of the ligand is :{pose_pdb.pdb_info().chain(i)}")
            for atm in range(1, pose_pdb.residue(i).natoms() + 1):
                x = round(pose_pdb.residue(i).xyz(atm)[0], 3)
                y = round(pose_pdb.residue(i).xyz(atm)[1], 3)
                z = round(pose_pdb.residue(i).xyz(atm)[2], 3)
                atom_name = pose_pdb.residue(i).atom_name(atm)
                # ligand_name = pose_pdb.residue(i).name()[-3:]
                ligand_name = pose_pdb.residue(i).name()
                # print(ligand_name)
                ligand_name = ligand_name.split(":")[0]
                # print(ligand_name)
                ligand_name = ligand_name.split("_")
                # print(ligand_name)
                if len(ligand_name) < 2:
                    ligand_name = ligand_name[0]
                else:
                    ligand_name = ligand_name[1]
                    # print(ligand_name) if ligand_name == "FUC": print("FUC FOUND >> NOT
                    # SAVED>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>") continue
                ligand_name = ligand_name.replace("_", "")
                ligand_name = ligand_name.strip()
                chain_name = pose_pdb.pdb_info().chain(i)
                num_chain = pose_pdb.pdb_info().pose2pdb(i)
                res_num = int(num_chain.split(" ")[0])
                save_pdb.save(file_path=save_path, x_cord=x, y_cord=y, z_cord=z, atom_name=atom_name,
                              ligand_name=ligand_name, chain_id=chain_name, residue_number=res_num, count=count)
                count += 1
            save_pdb.save_TER(file_path=save_path)
        save_pdb.save_end(file_path=save_path)


def check_ligands(pdb_template):
    # pose.residue(100).is_ligand
    # dump_pdb(pose,"/path/to/output_file.pdb")
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Checking Ligands")
    pose_pdb = pose_from_pdb(pdb_template)
    print(f"Total number of residue found for the pdb protein is : {pose_pdb.total_residue()}")
    for i in range(1, pose_pdb.total_residue() + 1):
        if pose_pdb.residue(i).is_ligand():
            # if pose_pdb.residue(i).name()[-3:] == casp_ligand:
            return True
