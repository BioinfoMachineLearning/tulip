"""
author: nabin 
timestamp: Thu Apr 21 2022 12:02 AM


** adapted from Ashwin's script **
"""
import os
import os.path as osp
import __main__

__main__.pymol_argv = ['pymol', '-qc']

import pymol

pymol.finish_launching()


def generate_complex(predicted, equibind_protein_path, equibind_ligand_path, save_path):
    os.makedirs(save_path, exist_ok=True)
    for pred in predicted:
        pred_idx = predicted.index(pred)
        pred = pred.split(".")[0]
        pdb_file = osp.join(equibind_protein_path, pred, "predicted_protein.pdb")
        ligand_file = osp.join(equibind_ligand_path, pred, 'lig_equibind_corrected.sdf')
        target_name = f"equibind_complex_pred_{pred_idx}.pdb"
        pdb_name = os.path.basename(pdb_file)
        ligand_name = os.path.basename(ligand_file)

        pymol.cmd.load(pdb_file, pdb_name)
        pymol.cmd.load(ligand_file, ligand_name)
        pymol.cmd.disable("all")
        pymol.cmd.enable(pdb_name)
        pymol.cmd.enable(ligand_name)
        print(pymol.cmd.get_names())
        print(target_name)
        save = osp.join(save_path, target_name)
        pymol.cmd.save(save)
        pymol.cmd.reinitialize()
    pymol.cmd.quit()



