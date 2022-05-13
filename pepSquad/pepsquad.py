import Bio.PDB as bio
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from operator import itemgetter

import numpy as np

import sys
import glob

from tqdm import tqdm

import itertools as it
source = "test"
source = sys.argv[1]
file_list = glob.glob(source+"/*.pdb")

l = 3
r = 2
b = 2.5
d = 0

l = int(sys.argv[2])
r = int(sys.argv[3])
b = float(sys.argv[4])
d = int(sys.argv[5])

def load_protein(source, name):
    parser = PDBParser()
    pdb = parser.get_structure(name,source)
    return pdb

def extract_lmer(proteins,i,l):
    residues = proteins.get_residues()
    list = []
    for res in residues:
        list.append(res.get_unpacked_list())
    return list[i:i+l]

def write_lmer_xyz(name, proteins,i,l):
    file = open(name+'.xyz', 'w')
    residues = proteins.get_residues()
    xyz = []
    list = []
    names = []
    
    for res in residues:
        list.append(res.get_unpacked_list())
        names.append(res.get_resname())

    for x in names[i:i+l]: print (x)
    print
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

#print atom.get_id(),list(atom.get_vector())
def extract_content(proteins):
    P = []
    for protein in proteins:
        atoms = protein.get_atoms()
        P_i = []
        for ca in atoms:
            P_i.append(ca.get_vector())
        P.append(P_i)
    return P

def extract_CAs(proteins):
    P = []
    for protein in proteins:
        atoms = protein.get_atoms()
        P_i = []
        for ca in atoms:
            if ca.get_id() == 'CA':
                P_i.append(ca.get_vector())
        P.append(P_i)
    return P

def get_id(file):
    #pdb = file[:-7]
    pdb = file[-8:-4]
    return pdb

def center_all(P):
    TP = []
    for p in P:
        p = np.array(p)
        center =  np.sum(p,axis=0) / float(len(p))
        p -= center
        TP.append(p)
    return TP

def get_all(l,r_sample, TP): return list(it.product(*[range(y-l+1) for y in [len(x) for x in [TP[i] for i in r_sample]]]))

def get_samples(l,r,TP):
    n = len(TP)
    samples = []
    for sample in  it.combinations(range(n), r):
        for item in get_all(l, sample, TP):
            lmers = []
            for i in range(len(sample)):
                lmers.append((sample[i],item[i]))
            samples.append(lmers)
    return samples

def get_lmers(l,sample,TP):
    lmers = []
    for (i,j) in sample:
        lmers.append(TP[i][j:j+l])
    return lmers

def convert(lvs):
    mvs = []
    for vector in lvs:
        x =[]
        for i in vector:
            x.append(i)
        mvs.append(np.array(x))
    return np.array(mvs)

def superimpose_pair(mol1, mol2):
    sel1 = np.array(mol1)
    sel2 = np.array(mol2)
    # check for consistency
    assert len(sel1) == len(sel2)
    L = len(sel1)
    assert L > 0

    # must always center the two proteins to avoid
    # affine transformations.  Center the two proteins
    # to their selections.

    COM1 = np.sum(sel1,axis=0) / float(L)
    COM2 = np.sum(sel2,axis=0) / float(L)
    sel1 -= COM1
    sel2 -= COM2
    csel1 = convert(sel1)
    csel2 = convert(sel2)

    # Initial residual, see Kabsch.
    E0 = np.sum( np.sum(csel1 * csel1,axis=0),axis=0) + np.sum( np.sum(csel2 * csel2,axis=0),axis=0)
    # This beautiful step provides the answer. V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!

    V, S, Wt = np.linalg.svd( np.dot( np.transpose(csel2), csel1))

    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    #U is simply V*Wt
    U = np.dot(V, Wt)

    # rotate and (not yet translate) the molecule
    sel2 -=COM2
    sel2 = np.dot(convert(sel2), U)

    # center the molecule
    sel1 = convert(sel1 - COM1)

    return (RMSD,sel1, sel2, U)

def pairwise_score(mols): # pairwise alignment
    r = len(mols)
    TRMSD = 0.0
    D = []
    for i in range(r):
        for j in range(i+1,r):
            (RMSD, cmol1, cmol2, U) = superimpose_pair(mols[i], mols[j])
            D.append([(i,j),RMSD,cmol1,cmol2, U])
            TRMSD += RMSD
    return TRMSD

def superimpose_samples_2(mols): # pairwise alignment
    r = len(mols)
    TRMSD = 0.0
    D = []
    aligned = []
    cmol = np.array([[0.0 for i in range(3)] for j in range(len(mols[0]))])
    #get pairwise rmsd
    for i in range(r):
        for j in range(i+1,r):
            (RMSD, mol1, mol2, U) = superimpose_pair(mols[i], mols[j])
            D.append([[i,j],RMSD,mol1,mol2, U])
            TRMSD += RMSD
    D_sorted = sorted(D, key=itemgetter(1))

    A = set()
    B = set()
    for d in D_sorted:
        A.update(d[0])
        if len(A-B)<1: continue
        elif len(A-B)==1:
            #print d[0], d[1], A-B, "Align ", list(A-B)[0], "to ", list(set(d[0]) - set(A-B))[0]
            j = list(A-B)[0]
            i = list(set(d[0]) - set(A-B))[0]
            (rmsd, mol1, mol2, rot) = superimpose_pair(mols[i],mols[j])
            cmol += mol2
            aligned.append(mol2)
        else:
            #print d[0], d[1], A-B, "Align ", d[0][0]," and ", d[0][1]
            i = d[0][0]
            j = d[0][1]
            (rmsd, mol1, mol2, rot) = superimpose_pair(mols[i],mols[j])
            cmol += mol1
            cmol += mol2
            aligned.append([mol1, mol2])
        TRMSD += rmsd
       
        B = set(A)
        if len(A)==r: break
    
    return (TRMSD, cmol/float(r))

def superimpose_samples_minstar(mols): # min star alignment: head: all
    r = len(mols)
    min_TRMSD = 999.0
    min_cmol = None
    for i in range(r):
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

def get_min_aln(consensus,tp):
    l = len(consensus)
    n = len(tp)
    min_rmsd = 999.0
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

def get_ball_size(mol):
    x = max([dim[0] for dim in mol])- min([dim[0] for dim in mol])
    y = max([dim[1] for dim in mol])- min([dim[1] for dim in mol])
    z = max([dim[2] for dim in mol])- min([dim[2] for dim in mol])
    return max([x,y,z])/float(2)

def write_xyz(name, mol): # for CA only
    file = open(name+".xyz", 'w')
    file.write(str(len(mol))+"\n\n")
    for x in mol:
        file.write("CA "+str(np.array(x))[1:-1]+"\n")

def get_feasible(consensus,TP):
    mols = []
    TRMSD = 0.0
    occ = []
    for tp in TP:
        (i, rmsd, mol) = get_min_aln(consensus, tp)
        occ.append(i)
        #print TP.index(tp), "MIN RMSD:", rmsd
        #print TP.index(tp),"th MOL"
        #for x in mol: #print x
        mols.append(mol)
        TRMSD += rmsd
    return (occ,TRMSD,mols)

def ptas(l,r,b,d,P):
    n  = len(P)
    min_trmsd = 999.0
    min_feasible = None
    min_consensus = None
    min_occ = None
    ##############
    #   Steps
    #       1. Coincide centers of P
    #       2. For every r l-length motif
    #           3. Get list of vectors from sample
    #           4. Get minimum alignment of  mol_sample
    #           5.

    #1. Coincide centers of P
    #TP = coincide_centers(P)
    TP = center_all(P)

    #2. For every r l-length motif
    samples = get_samples(l,r,TP)
    m = len(samples)
    
    #print m
    count = 0
    for sample in tqdm(samples, desc="Calculating..."):
        #3. Get list of vectors from sample
        mol_sample = get_lmers(l,sample,TP)
        ball_size = max([get_ball_size(mol) for mol in mol_sample])
        if ball_size< b: #opt: pruning before alignment in max ballsize outputs consensus less than obtd ball size
            count +=1
            #4. Get minimum alignment of  mol_sample
            (rmsd,consensus) = superimpose_samples(mol_sample)
            (occ,trmsd,feasible) = get_feasible(consensus,TP)
            #print "TRMSD:", trmsd, occ
            if min_trmsd > trmsd:
            #if min_trmsd > pairwise_score(feasible):
                min_trmsd = trmsd
                #min_trmsd = pairwise_score(feasible)
                min_consensus = consensus
                min_feasible = feasible
                min_occ = occ
                min_sample = sample

    print ("EVALUATED: ", count, '/', m, '=', (float(count)/m) * 100.0)

    if count > 0:
        print ("PAIRWISE RMSD:", pairwise_score(min_feasible))
        return (min_occ, min_trmsd, min_consensus, min_feasible, min_sample)
    else:
        print ("No Structural Motifs Found!")
        sys.exit()

def center(mol):
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

def align_back_to_ca_s(filename,CA_s):
    P = []
    P_CAs = []
    names =[]
    ns = []
    fn = []
    for file in glob.glob(filename[:11]+'*'+filename[12:]+".xyz"):
        p =[]
        protein = load_xyz(file)
        fn.append(file)
        
        ns.append(protein[0])
        names.append(protein[1])
        for atom in protein[2]:
            p.append(np.array([float(x) for x in atom]))
        P.append(p)

    for p_names,vectors in zip(names,P):
        P_CAs.append([v for (n,v) in zip(p_names,vectors) if n=='CA'])

    for i in range(len(P_CAs)):
        (d, mol1, mol2, rot) =  superimpose_pair(CA_s[i],P_CAs[i])
        
        aligned = np.dot(center(P[i]),rot)
        file = open("ALN_"+fn[i],'w')
    
        file.write(str(ns[i])+'\n')
        for j in range(ns[i]):
            file.write(names[i][j]+' '+str(aligned[j])[1:-1]+'\n')
        file.close()

if __name__ == "__main__":
    proteins = []
    pdb_ids = []

    for file in file_list:
        id = get_id(file)
        pdb_ids.append(id)
        proteins.append(load_protein(file,id))

    backbones = extract_content(proteins)
    #print len(backbones)

    #for p in backbones: print len(p),
    #   print
    P = extract_CAs(proteins)
    #print len(P)
    #for p in P: print len(p),
    (occ, rmsd, pattern, occurence, sample) = ptas(l,r,b,d,P)
    print ("SAMPLE:", sample)
    print ("BEST_RMSD:", rmsd)
    print ("BEST_OCCURENCE:", occ)
    print ("BALL SIZE:", get_ball_size(pattern))

    family = source[source.rfind('/')+1:]
    name = family+'_'+str(l)+'_'+str(r)+'_'+str(b).replace('.','-')+'_'+str(d)
    print ("CONSENSUS:")
    write_xyz("CONSENSUS_"+name,center(pattern))
    for x in pattern:
        print ("CA", x)
    print ("OCCURRENCE:")
    for i in range(len(occ)):
        print (i+1,'\t',pdb_ids[i],'\t',occ[i]+1,'\t',occ[i]+l,'\t',)
        write_lmer_xyz("OCCUR_BACK_"+str(i)+'_'+name,proteins[i],occ[i], l)
        write_xyz("OCCUR_CA_"+str(i)+name, center(occurence[i]))

        align_back_to_ca_s("OCCUR_BACK_"+str(i)+'_'+name, occurence)
