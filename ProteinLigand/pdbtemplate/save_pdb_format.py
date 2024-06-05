"""
author: nabin 
timestamp: Wed Apr 20 2022 3:56 PM
"""
AUTHOR_NAME = "BML"

def save_titles(file_path, target, remark, method, model_number):
    with open(file_path, 'a') as fi:
        fi.write('PFRMAT')
        fi.write('  ')
        fi.write('TS')
        fi.write('\n')
        fi.write('TARGET')
        fi.write('  ')
        fi.write(str(target))
        fi.write('\n')
        fi.write('AUTHOR')
        fi.write('  ')
        fi.write(AUTHOR_NAME)
        fi.write('\n')
        fi.write('REMARK')
        fi.write('  ')
        fi.write(remark)
        fi.write('\n')
        fi.write('METHOD')
        fi.write('  ')
        fi.write(method)
        fi.write('\n')
        fi.write('MODEL')
        fi.write('  ')
        fi.write(str(model_number))
        fi.write('\n')
        fi.write('PARENT')
        fi.write('  ')
        fi.write(' ')
        fi.write('\n')
    fi.close()


def save(file_path, x_cord, y_cord, z_cord, atom_name, ligand_name, chain_id, residue_number, count):
    with open(file_path, 'a') as fi:
        fi.write('HETATM')
        fi.write('  ')
        fi.write(str(count).rjust(3))
        fi.write('  ')
        fi.write(atom_name.ljust(4))
        fi.write(ligand_name.rjust(3))
        fi.write(' ')
        fi.write(chain_id)
        fi.write(str(residue_number).rjust(4))
        fi.write('    ')
        fi.write(str(x_cord).rjust(8))
        fi.write(str(y_cord).rjust(8))
        fi.write(str(z_cord).rjust(8))
        fi.write(str(1.00).rjust(5))
        fi.write(str(0.00).rjust(5))
        fi.write('           ')
        fi.write(atom_name[0:1].rjust(1))
        fi.write('  ')
        fi.write('\n')


def save_TER(file_path):
    with open(file_path, 'a') as fi:
        fi.write('TER')
        fi.write('\n')


def save_end(file_path):
    with open(file_path, 'a') as fi:
        fi.write('END')
        fi.write('\n')
