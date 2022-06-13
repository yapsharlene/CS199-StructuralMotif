from pdb_extract import *
from tqdm import tqdm
import sys
import itertools as it
import numpy as np
import multiprocessing
import time

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

# from pepsquad
def convert(lvs):
    mvs = []
    for vector in lvs:
        x =[]
        for i in vector:
            x.append(i)
        mvs.append(np.array(x))
    return np.array(mvs)

# from pepsquad
def superimpose_pair(mol1, mol2):
    sel1 = np.array(mol1)
    sel2 = np.array(mol2)
    
    assert len(sel1) == len(sel2)
    L = len(sel1)
    assert L > 0

    COM1 = np.sum(sel1,axis=0) / float(L)
    COM2 = np.sum(sel2,axis=0) / float(L)
    sel1 -= COM1
    sel2 -= COM2
    csel1 = convert(sel1)
    csel2 = convert(sel2)

    E0 = np.sum( np.sum(csel1 * csel1,axis=0),axis=0) + np.sum( np.sum(csel2 * csel2,axis=0),axis=0)

    V, S, Wt = np.linalg.svd( np.dot( np.transpose(csel2), csel1))

    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    U = np.dot(V, Wt)

    sel2 -=COM2
    sel2 = np.dot(convert(sel2), U)

    sel1 = convert(sel1 - COM1)

    return (RMSD,sel1, sel2, U)

# from pepsquad
def pairwise_score(mols): # pairwise alignment
    r = len(mols)
    TRMSD = 0.0
    D = []
    for i in tqdm(range(r), desc = "getting pairwise score"):
        for j in range(i+1,r):
            (RMSD, cmol1, cmol2, U) = superimpose_pair(mols[i], mols[j])
            D.append([(i,j),RMSD,cmol1,cmol2, U])
            TRMSD += RMSD
    return TRMSD

# from pepsquad
def superimpose_samples(mols): # min star alignment: head: all
    r = len(mols)
    min_TRMSD = 999.0
    min_cmol = None
    for i in range(1):
        TRMSD = 0.0
        cmol = convert(mols[i])
        for mol2 in mols[:i]+mols[i+1:]:
            (RMSD, cmol1, cmol2, U) = superimpose_pair(mols[i], mol2)
            TRMSD += RMSD
            cmol += cmol2
        if min_TRMSD > TRMSD:
            min_TRMSD = TRMSD
            min_cmol = cmol
    return (min_TRMSD, min_cmol/float(r))

# from pepsquad
def get_min_aln(consensus,tp):
    l = len(consensus)
    n = len(tp)
    min_rmsd = sys.maxsize
    cmol = None
    j = 0
    for mol in [tp[i:i+l] for i in range(n-l)]:
        j +=1
        (rmsd, mol1, mol2, U) = superimpose_pair(consensus, mol)
        if min_rmsd > rmsd:
            min_rmsd = rmsd
            min_j = j
            cmol = mol2
    return (min_j, min_rmsd,cmol)

# from pepsquad
def get_feasible(consensus,TP):
    mols = []
    TRMSD = 0.0
    occ = []
    for tp in TP:
        (i, rmsd, mol) = get_min_aln(consensus, tp)
        occ.append(i)
        mols.append(mol)
        TRMSD += rmsd
    return (occ,TRMSD,mols)

# superimposition process for parallelization
def impose_step(data):
    sample_ranges = data[0]
    centered_structs = data[1]
    min_trmsd = data[2]
    length = data[3]
    b = data[4]
    ID = data[5]

    count = 0
    min_feasible = None
    min_consensus = None
    min_occ = None
    min_ball_size = None
    min_sample = None
    
    for sample_range in tqdm(sample_ranges, desc = "superimposing ID: " + ID):
        # extract the l-length motifs
        motifs = extract(centered_structs, sample_range, length)

        # calculate ball size 
        ball_size = get_ball_size(motifs)

        if ball_size < b:
            count += 1

            (rmsd,consensus) = superimpose_samples(motifs)
            (occ,trmsd,feasible) = get_feasible(consensus,centered_structs)

            if min_trmsd > trmsd:
                min_trmsd = trmsd
                min_consensus = consensus
                min_feasible = feasible
                min_occ = occ
                min_sample = sample_range
                min_ball_size = ball_size

    return [count,min_trmsd,min_consensus,min_feasible,min_occ,min_sample,min_ball_size]

# === MAIN === #
if __name__ == "__main__":
    
    tic = time.time()
    
    BENCHMARK_LENGTH = 3
    r = 2
    b = 2.5
    d = 0
    min_trmsd = sys.maxsize

    # Get all protein structs
    structs = get_structures("/datasets","/pep_conantokin/")

    # (1) Fix P_1, translate other proteins to make centroids coincide
    fixed_struct = structs[0]

    # centering proteins
    centered_structs = center(structs)

    # getting samples (all r l-length motif)
    sample_ranges = sample(centered_structs, r, BENCHMARK_LENGTH)

    print(len(sample_ranges))

    # parallelization step
    num_processors = 1
    div = int(len(sample_ranges)/num_processors)
    pool = multiprocessing.Pool(num_processors)
    results = pool.map(impose_step,[[sample_ranges[i:i+div],centered_structs,min_trmsd,BENCHMARK_LENGTH,b,str(int(i/div))] for i in range(0,len(sample_ranges),div)])
    pool.close()

    print(len(results))
    
    # get total count and result w/ smallest min_trmsd (first two elements)
    minimum = min(results, key = lambda x: x[1])
    tot_count = sum(i[0] for i in results)
    minimum[0] = tot_count
    
    print("EVALUATED: ", tot_count, '/', len(sample_ranges), '=', (float(tot_count)/len(sample_ranges)) * 100.0)

    # min_feasible is at index 3
    if tot_count > 0:
        print("PAIRWISE RMSD:", pairwise_score(minimum[3]))
    else:
        print("No Structural Motifs Found!")
    
    # [count,min_trmsd,min_consensus,min_feasible,min_occ,min_sample,min_ball_size]
    print("SAMPLE: ", minimum[5])
    print("BEST_RMSD: ", minimum[1])
    print("BEST_OCCURENCE", minimum[4])
    print("BALL SIZE, ", minimum[6])

    toc = time.time()
    print("Program done in {:.4f} seconds".format(toc-tic))