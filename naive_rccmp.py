from pdb_input import *

# === NAIVE RCCMP === #
# INPUT:
# - n protein structures (P_1, P_2, ..., P_n)
# - length of RCCM l
# - C number of residues
# - r samples


# === MAIN === #

structs = get_structures("PDB files\Ref2")
for struct in structs:
    grab_ca_coords(struct)






