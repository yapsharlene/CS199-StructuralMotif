from test import *
import os
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
def get_range(proteins, r_sample, l):
    # ranges = []
    # for i in r_sample:
    #     ranges.append(range(len(proteins[i])-l+1))
    # return ranges
    return list(it.product(*[range(y-l+1) for y in [len(x) for x in [proteins[i] for i in r_sample]]]))

# gets the ranges of the motifs to be sampled
def sample(proteins, r, l):
    samples = []    
    for indices in it.combinations(range(len(proteins)),r):
        for ranges in get_range(proteins, indices, l):
            sample = []
            for i in range(len(indices)):
                sample.append((indices[i],ranges[i]))
            samples.append(sample)
    return samples

# extracts the l-length motif given range
def extract(proteins, idx, l):
    motif = []
    for (i,j) in idx:
        motif.append(proteins[i][j:j+l])
    return motif

# calculates the largest ball size possible given the motifs
def get_ball_size(motifs):
    ball_sizes = []
    for motif in motifs:
        x = max([dim[0] for dim in motif])- min([dim[0] for dim in motif])
        y = max([dim[1] for dim in motif])- min([dim[1] for dim in motif])
        z = max([dim[2] for dim in motif])- min([dim[2] for dim in motif])
        ball_sizes.append(max([x,y,z])/float(2))   
    return max(ball_sizes)
# === MAIN === #

BENCHMARK_LENGTH = 3
r = 2

# Get all protein structs
# structs = get_structures("PDB files\Ref1")
structs = get_structures("/datasets","/pdb/")

# for i in structs[1]:
#     print(i)

# (1) Fix P_1, translate other proteins to make centroids coincide
fixed_struct = structs[0]

# centering proteins
centered_structs = center(structs)
# for i in range(1,len(structs)):
#     centroid = get_centroid(structs[i])
#     for j in range(len(structs[i])):
#         for k in range(3):
#             structs[i][j][k] -= centroid[k]

# print("======CENTERED!======")
# for i in centered_structs[1]:
#     print(i)

# getting samples (all r l-length motif)
sample_ranges = sample(centered_structs, r, BENCHMARK_LENGTH)

for sample_range in sample_ranges:
    # extract the l-length motifs
    motifs = extract(centered_structs, sample_range, BENCHMARK_LENGTH)

    # calculate ball size 
    ball_size = get_ball_size(motifs)

#     superimposer = Bio.PDB.Superimposer()
#     for i in range(1, len(structs)):
#         # print_ca_coords(get_ca_atoms(structs[i], 3))
#         superimposer.set_atoms(get_ca_atoms_of_range(fixed_struct, BENCHMARK_LENGTH), get_ca_atoms_of_range(structs[i], BENCHMARK_LENGTH))
#         superimposer.apply(get_ca_atoms_of_range(structs[i], BENCHMARK_LENGTH))
#         # print(superimposer.rms)
#         # print(superimposer.rotran)
#         print("Superimposed alpha carbons of {0} to be as close as possible to {1}".format(structs[i].id, fixed_struct.id))

# (2) Select a length-l compact motif u1, u2, …, ur
# where ui is a motif of some Pj and x is the total number of samples
        

# (3) Select a transformation 𝜏2𝜏3,...,𝜏r

# (3a) Find the median for each discrete rigid transformation u

# (3b) Find the compact motif that minimizes the cMSD distance vi where i = 1, 2, … n

# (3c) Compute the objective function value c(u)
