from pdb_input import *

import itertools as it

# INPUT:
# - n protein structs (P_1, P_2, ..., P_n)
# - length of RCCM l
# - C number of residues
# - r samples

# coincides proteins with their centroid
def center(proteins):
    centered = []
    for protein in proteins:
        centroid = get_centroid(protein)
        for j in range(len(protein)):
            for k in range(3):
                protein[j][k] -= centroid[k]
        centered.append(protein)
    return centered

# gets valid range in extracting l-mers 
def get_range(proteins, r, l):
    ranges = []
    for i in range(r):
        ranges.append(range(len(proteins[i])-l+1))
    return ranges

# gets all r l-length motifs
def sample(proteins, r, l):
    samples = []    
    for indices in it.combinations(range(len(proteins),r)):
        for range in get_range(proteins, indices, l):
            sample = []
            for i in range(len(indices)):
                sample.append((indices[i],range[i]))
            samples.append(sample)
    return sample

# === MAIN === #

BENCHMARK_LENGTH = 3

# Get all protein structs
structs = get_structures("PDB files\Ref1")

# (1) Fix P_1, translate other proteins to make centroids coincide
fixed_struct = structs[0]

# centering proteins
centered_structs = center(structs[1:])
# for i in range(1,len(structs)):
#     centroid = get_centroid(structs[i])
#     for j in range(len(structs[i])):
#         for k in range(3):
#             structs[i][j][k] -= centroid[k]

# getting samples
sample_residues = sample(centered_structs)

superimposer = Bio.PDB.Superimposer()
for i in range(1, len(structs)):
    # print_ca_coords(get_ca_atoms(structs[i], 3))
    superimposer.set_atoms(get_ca_atoms_of_range(fixed_struct, BENCHMARK_LENGTH), get_ca_atoms_of_range(structs[i], BENCHMARK_LENGTH))
    superimposer.apply(get_ca_atoms_of_range(structs[i], BENCHMARK_LENGTH))
    # print(superimposer.rms)
    # print(superimposer.rotran)
    print("Superimposed alpha carbons of {0} to be as close as possible to {1}".format(structs[i].id, fixed_struct.id))

# (2) Select a length-l compact motif u1, u2, …, ur
# where ui is a motif of some Pj and x is the total number of samples
        

# (3) Select a transformation 𝜏2𝜏3,...,𝜏r

# (3a) Find the median for each discrete rigid transformation u

# (3b) Find the compact motif that minimizes the cMSD distance vi where i = 1, 2, … n

# (3c) Compute the objective function value c(u)
