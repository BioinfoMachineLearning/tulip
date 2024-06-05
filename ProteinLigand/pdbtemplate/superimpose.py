"""
author: nabin 
timestamp: Wed Apr 20 2022 4:26 PM

1. runs chimera's matchmaker in no GUI mode to superimpose predicted structure and PDB template structure that has provided ligands in it

2. saves the superimposed structure relative to PDB template



"""
import os
import subprocess


def superimpose_structures(chimera_path, pdb_template, predicted_structure, save_path, idx):
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Superimposing/Overlay structures")
    os.makedirs(save_path, exist_ok=True)
    super_impose_cmd_path = f"{save_path}/superimpose_struc.cmd"
    chimera_scripts = open(super_impose_cmd_path, 'w')
    chimera_scripts.write('open ' + predicted_structure + '\n'
                                                          'open ' + pdb_template + '\n'
                                                                                   'mm ' + '#0 ' + '#1' '\n'
                                                                                                   'select ' + '#0,1' '\n'
                                                                                                               'write ' + 'format ' + 'pdb ' + 'selected ' + 'relative ' + '0 ' + '1 ' + f'{save_path}/matchmaker_saved_model_{idx}.pdb')
    chimera_scripts.close()
    command_complete = False
    while not command_complete:
        try:
            cmd = [chimera_path, '--nogui', chimera_scripts.name]
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            stdout = result.stdout
            stderr = result.stderr
            return_code = result.returncode
            if return_code == 0:
                print("Superimpose executed successfully.")
                # print(stdout)
            else:
                print(f"Superimpose failed with exit code {return_code}.")
                print("Standard Error:")
                print(stderr)
            command_complete = True
        except FileNotFoundError as error:
            raise error
    
    # os.remove(chimera_scripts.name)
    if os.path.exists(chimera_scripts.name):
        os.remove(chimera_scripts.name)
        print(f"File '{chimera_scripts.name}' has been removed.")
    else:
        print(f"File '{chimera_scripts.name}' not found.")
