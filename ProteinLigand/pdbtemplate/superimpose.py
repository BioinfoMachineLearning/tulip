"""
author: nabin 
timestamp: Wed Apr 20 2022 4:26 PM

1. runs chimera's matchmaker in no GUI mode to superimpose predicted structure and PDB template structure which has
CASP15 provided ligands in it

2. saves the superimposed structure relative to PDB template


model is saved by name matchmaker_saved_model.pdb

"""
import os
import subprocess


def superimpose_structures(chimera_path, pdb_template, predicted_structure, save_path, idx):
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Superimposing/Overlay structures")
    os.makedirs(save_path, exist_ok=True)
    chimera_scripts = open('superimpose_struc.cmd', 'w')
    chimera_scripts.write('open ' + predicted_structure + '\n'
                                                          'open ' + pdb_template + '\n'
                                                                                   'mm ' + '#0 ' + '#1' '\n'
                                                                                                   'select ' + '#0,1' '\n'
                                                                                                               'write ' + 'format ' + 'pdb ' + 'selected ' + 'relative ' + '0 ' + '1 ' + f'{save_path}/matchmaker_saved_model_{idx}.pdb')
    chimera_scripts.close()
    command_complete = False
    while not command_complete:
        try:
            subprocess.run([chimera_path, '--nogui', chimera_scripts.name])
            command_complete = True
        except FileNotFoundError as error:
            raise error

    os.remove(chimera_scripts.name)
