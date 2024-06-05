
import os

def sep_identical_ligs(PDB_template_ligand, idx):

    ligand_names = list()

    with open(PDB_template_ligand, 'r') as file:
        contents = file.readlines()

    current_section = list()
    count = 0
    lig_parent_dir = os.path.dirname(PDB_template_ligand)
    lig_name = os.path.basename(PDB_template_ligand)
    lig_name = lig_name.split(".")[0]
    lig_name = lig_name.split("_")
    lig_idx = lig_name[1]
    lig_name = lig_name[0]
    for line in contents:
        current_section.append(line)

        # Check if the line contains "TER"
        if "TER\n" in line:
            with open(f'{lig_parent_dir}/{lig_name}{count}_{lig_idx}.pdb', 'w') as output_file:
                output_file.writelines(current_section)
            current_section = list() 
            ligand_names.append(f'{lig_name}{count}')
            count += 1
        if "END\n" in line:
            break
    
    return ligand_names