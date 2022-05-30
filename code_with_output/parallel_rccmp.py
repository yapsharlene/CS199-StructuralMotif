from test import *
from tqdm import tqdm
import sys
import itertools as it
import numpy as np
import multiprocessing
import time
import glob

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

def par_superimpose_samples(data):
    mols = data[0]
    i = data[1]
    min_TRMSD = sys.maxsize
    min_cmol = None
    TRMSD = 0.0
    cmol = convert(mols[i])
    for mol2 in mols[:i]+mols[i+1:]:
        (RMSD, cmol1, cmol2, U) = superimpose_pair(mols[i], mol2)
        TRMSD += RMSD
        cmol += cmol2
    if min_TRMSD > TRMSD:
        min_TRMSD = TRMSD
        min_cmol = cmol
    return (min_TRMSD, min_cmol)

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

    # processors = 2
    # pool = multiprocessing.Pool(processors)
    # results = pool.map(par_superimpose_samples,[[mols,i] for i in range(1)])
    # minimum = min(results, key = lambda x: x[0])
    # return (minimum[0], minimum[1]/float(r))

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

def center_2(mol):
    return mol - np.sum(mol,axis=0)/float(len(mol))

def load_xyz(filename):
    xyzs = []
    names = []
    file = open(filename, 'r')
    n = int(file.readline())
    for line in file.readlines():
        xyzs.append(line.split(' ')[1:-1])
        names.append(line.split(' ')[0])
    file.close()
    return (n, names, xyzs)

def write_xyz(name, mol): # for CA only
    file = open(name+".xyz", 'w')
    file.write(str(len(mol))+"\n\n")
    for x in mol:
        file.write("CA "+str(np.array(x))[1:-1]+"\n")

def align_back_to_ca_s(filename,CA_s):
    P = []
    P_CAs = []
    names =[]
    ns = []
    fn = []
    #print filename[:11]+'*'+filename[12:]+".xyz"
    #print len(glob.glob(filename[:11]+'*'+filename[12:]))
    for file in glob.glob(filename[:11]+'*'+filename[12:]+".xyz"):
        p =[]
        protein = load_xyz(file)
        fn.append(file)
        #print file
        
        ns.append(protein[0])
        names.append(protein[1])
        for atom in protein[2]:
            p.append(np.array([float(x) for x in atom]))
        P.append(p)

    for p_names,vectors in zip(names,P):
        P_CAs.append([v for (n,v) in zip(p_names,vectors) if n=='CA'])

    for i in range(len(P_CAs)):
        (d, mol1, mol2, rot) =  superimpose_pair(CA_s[i],P_CAs[i])
        #print d
        
        aligned = np.dot(center_2(P[i]),rot)
        file = open("ALN_"+fn[i],'w')
        #print "ALN_"+fn[i]
    
        file.write(str(ns[i])+'\n')
        for j in range(ns[i]):
            file.write(names[i][j]+' '+str(aligned[j])[1:-1]+'\n')
        file.close()

def get_id(file):
    #pdb = file[:-7]
    pdb = file[-8:-4]
    return pdb

def load_protein(source, name):
    parser = PDBParser()
    pdb = parser.get_structure(name,source)
    return pdb

def write_lmer_xyz(name, proteins,i,l):
    file = open(name+'.xyz', 'w')
    residues = proteins.get_residues()
    xyz = []
    list = []
    names = []

    for res in residues:
        list.append(res.get_unpacked_list())
        names.append(res.get_resname())

    for x in names[i:i+l]: print(x)
    count = 0
    line = ""
    for res in list[i:i+l]:
        count += len(res)
        for atom in res:
            line = line+str(atom.get_id())+' '
            for i in atom.get_vector(): line = line+str(i)+' '
            line  = line+'\n'

    file.write(str(count)+'\n')
    file.write(line)

# === MAIN === #
if __name__ == "__main__":
    tic = time.time()
    BENCHMARK_LENGTH = 3
    r = 2
    b = 2.5
    d = 0
    min_trmsd = sys.maxsize

    # Get all protein structs
    # structs = get_structures("PDB files\Ref1")
    structs = get_structures("/datasets","/test/")

    # for struct in structs:
    #     print(len(struct))

    # (1) Fix P_1, translate other proteins to make centroids coincide
    fixed_struct = structs[0]

    # centering proteins
    centered_structs = center(structs)

    # getting samples (all r l-length motif)
    sample_ranges = sample(centered_structs, r, BENCHMARK_LENGTH)

    num_processors = 10
    div = int(len(sample_ranges)/num_processors)
    pool = multiprocessing.Pool(num_processors)
    results = pool.map(impose_step,[[sample_ranges[i:i+div],centered_structs,min_trmsd,BENCHMARK_LENGTH,b,str(int(i/div))] for i in range(0,len(sample_ranges),div)])
    pool.close()

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
    print("BEST_OCCURENCE: ", minimum[4])
    print("BALL SIZE: ", minimum[6])
    print("PATTERN: ", minimum[2])

    toc = time.time()

    print("Program done in {:.4f} seconds".format(toc-tic))

    source = "test"
    file_list = glob.glob("datasets/test/*.pdb")
    
    print(len(file_list))

    proteins = []
    pdb_ids = []

    for file in file_list:
        id = get_id(file)
        pdb_ids.append(id)
        proteins.append(load_protein(file,id))

    family = source[source.rfind('/')+1:]
    name = family+'_'+str(BENCHMARK_LENGTH)+'_'+str(r)+'_'+str(b).replace('.','-')+'_'+str(d)
    
    print("CONSENSUS:")
    write_xyz("CONSENSUS_"+name,center_2(minimum[2]))
    for x in minimum[2]:
        print("CA", x)

    print("OCCURENCE:")
    for i in range(len(minimum[4])):
        print(i+1,'\t',pdb_ids[i],'\t',minimum[4][i]+1,'\t',minimum[4][i]+BENCHMARK_LENGTH,'\t')
        write_lmer_xyz("OCCUR_BACK_"+str(i)+'_'+name,proteins[i],minimum[4][i], BENCHMARK_LENGTH)
        write_xyz("OCCUR_CA_"+str(i)+'_'+name, center_2(minimum[3][i]))
        align_back_to_ca_s("OCCUR_BACK_"+str(i)+'_'+name, minimum[3])
