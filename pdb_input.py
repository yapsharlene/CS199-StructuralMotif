import os
import warnings
from Bio.PDB.PDBParser import PDBParser

# === FUNCTIONS === #

# Get all structures in folder and return in list pdb_structs
def get_structures(folder):
    parser = PDBParser(PERMISSIVE=1)
    directory = folder

    pdb_structs = []
    id_num = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # remove general PDB warnings

        for file in os.listdir(directory):
            f = os.path.join(directory, file)
            if os.path.isfile(f):
                structure_id = os.path.splitext(file)[0]
                pdb_structs.append(parser.get_structure(structure_id, f))
                id_num += 1

    return pdb_structs

# Given a structure, print out the xyz-coordinates
# of each carbon atom
def print_ca_coords(struct):
    print('{:10}|{:10}|{:10}'.format("X","Y","Z"))
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_id() == "CA":
                        x,y,z = atom.get_coord()
                        print('{:<10.3f}|{:<10.3f}|{:<10.3f}'.format(x, y, z))
                        
# Given a structure, return a list of xyz-coordinates
# of each carbon atom
def grab_ca_coords(struct):
    atoms = []
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_id() == "CA":
                        atoms.append(atom.get_coord())
    print("Grabbed {0} coordinates from {1}".format(len(atoms), struct.id))
    return atoms