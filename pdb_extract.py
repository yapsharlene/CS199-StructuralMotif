import os
import itertools as it

# has PDBConstruction warnings; ignore by adding -W ignore in terminal
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)

# getting the filepath of the datasets
# ds_folder = "/datasets"
# refnum = "/pdb/"

def get_structures(ds_folder,refnum):
	path = os.getcwd() + ds_folder + refnum

	residues = []

	# process each pdb file
	for filename in os.listdir(path):
		structure = parser.get_structure(filename[:4], path + filename)
		coords = []

		# getting all CA coordinates of a structure
		for model in structure.get_list():
			for chain in model.get_list():
				for residue in chain.get_list():
					if residue.has_id("CA"):
						ca = residue["CA"]
						coords.append(ca.get_coord())

		residues.append(coords)
	return residues

# Return the centroid of a given atom list
def get_centroid(atom_list):
    centroid = [0, 0, 0]
    for atom in atom_list:
        centroid[0] += atom[0]
        centroid[1] += atom[1]
        centroid[2] += atom[2]

    for i in range(3):
        centroid[i] /= float(len(atom_list))

    return (centroid[0], centroid[1], centroid[2])