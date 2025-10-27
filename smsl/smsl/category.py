import numpy as np
import pandas as pd

_bond_data_from_charmm = {
"DG": f"""
BOND P    O1P       P    O2P       P     O5'
BOND O5'  C5'       C5'  C4'       C4'  O4'       C4'  C3'       O4'  C1'
BOND C1'  N9        C1'  C2'       N9   C4        N9   C8        C4   N3
BOND C2   N2        C2   N1        N2   H21
BOND N2   H22       N1   H1        N1   C6        C6   C5
BOND C5   N7        C2'  C3'       C3'  O3'       O3'  +P
BOND C2'  O2'       O2'  H2'
BOND C1'  H1'       C2'  H2''      C3'  H3'       C4'  H4'       C5'  H5'
BOND C5'  H5''      C8   H8
""",
"DA": f"""
BOND P    O1P       P    O2P       P     O5'
BOND O5'  C5'       C5'  C4'       C4'  O4'       C4'  C3'       O4'  C1'
BOND C1'  N9        C1'  C2'       N9   C4        N9   C8        C4   N3
BOND C2   N2        C2   N1        N2   H21
BOND N2   H22       N1   H1        N1   C6        C6   C5
BOND C5   N7        C2'  C3'       C3'  O3'       O3'  +P
BOND C2'  O2'       O2'  H2'
BOND C1'  H1'       C2'  H2''      C3'  H3'       C4'  H4'       C5'  H5'
BOND C5'  H5''      C8   H8
""",
"DC": f"""
BOND P    O1P       P    O2P       P     O5'
BOND O5'  C5'       C5'  C4'       C4'  O4'       C4'  C3'       O4'  C1'
BOND C1'  N1        C1'  C2'       N1   C2        N1   C6
BOND C2   N3        C4   N4        N4   H41       N4   H42
BOND C4   C5        C2'  C3'       C3'  O3'       O3'  +P
BOND C2'  O2'       O2'  H2'
BOND C1'  H1'       C2'  H2''      C3'  H3'       C4'  H4'       C5'  H5'
BOND C5'  H5''      C5   H5        C6   H6
""",
"DT": f"""
BOND P    O1P       P    O2P       P     O5'
BOND O5'  C5'       C5'  C4'       C4'  O4'       C4'  C3'       O4'  C1'
BOND C1'  N1        C1'  C2'       N1   C2        N1   C6
BOND C2   N3        C4   N4        N4   H41       N4   H42
BOND C4   C5        C2'  C3'       C3'  O3'       O3'  +P
BOND C2'  O2'       O2'  H2'
BOND C1'  H1'       C2'  H2''      C3'  H3'       C4'  H4'       C5'  H5'
BOND C5'  H5''      C5   H5        C6   H6
""",
}

def check_atomname_ij_is_moiety(df, moiety, times, type_na):
    is_moiety = get_is_moiety(df, moiety, type_na)
    n_isin = is_moiety.sum(axis=1)
    return n_isin==times

def get_is_moiety(df, moiety, type_na):
    moieties2atomname = {
        'R': "'",
        'P': "P",
    }
    atomname_ij = df[['Atomname_i','Atomname_j']].values
    strand_ij   = df[['Strand_i','Strand_j']].values
    if moiety=='B':
        is_moieties = []
        moieties = ['D', 'R', 'P'] if type_na=='dsdna' else ['M', 'D', 'R', 'P']
        for moiety in moieties:
            is_moiety = get_is_moiety(df, moiety, type_na)
            is_moieties.append(is_moiety)
        is_moiety = np.logical_not(np.logical_or.reduce(is_moieties))
    elif moiety=='M':
        strand_id = 'STRAND2'
        is_moiety = np.vectorize(lambda s: strand_id==s)(strand_ij)
    elif moiety=='D':
        strand_id = 'STRAND3'
        is_moiety = np.vectorize(lambda s: strand_id==s)(strand_ij)
    elif moiety=='R':
        atomname = "'"
        is_moiety = np.vectorize(lambda s: atomname in s)(atomname_ij)
    elif moiety=='P':
        atomname = "P"
        is_moiety = np.vectorize(lambda s: atomname in s)(atomname_ij)
    return is_moiety

def check_resid_ij_is_nucleobase(df, base_type, times, sequence):
    resid_ij = df[['Resid_i','Resid_j']].values
    is_moiety = np.vectorize(lambda resid: sequence[resid]==base_type)(resid_ij)
    n_isin = is_moiety.sum(axis=1)
    return n_isin==times

def GetResidLayer(resids, df_tetrad_geometry):
    layer_from_resid = np.empty_like(resids, dtype='str')
    for layer in df_tetrad_geometry.index:
        is_layer = np.where(np.isin(resids, df_tetrad_geometry.T[layer]),layer, '') # 'TOP', '', 'TOP', '' 
        layer_from_resid = np.char.add(layer_from_resid, is_layer)
    return layer_from_resid

def GetResidStrand(resids, df_tetrad_geometry):
    strand_from_resid = np.empty_like(resids, dtype='str')
    for strand in df_tetrad_geometry.columns:
        is_layer = np.where(np.isin(resids, df_tetrad_geometry[strand]),strand, '') # 'QIII', '', 'QII', '' 
        strand_from_resid = np.char.add(strand_from_resid, is_layer)
    return strand_from_resid

def check_resid_ij_same_layer(df, df_tetrad_geometry):
    layer_i = GetResidLayer(df['Resid_i'], df_tetrad_geometry)
    layer_j = GetResidLayer(df['Resid_j'], df_tetrad_geometry)
    return layer_i==layer_j

def check_resid_ij_same_strand(df, df_tetrad_geometry):
    strand_i = GetResidStrand(df['Resid_i'], df_tetrad_geometry)
    strand_j = GetResidStrand(df['Resid_j'], df_tetrad_geometry)
    return strand_i==strand_j

## Bonad And Angle
def GetBondAnglePair(bond_pairs=_bond_data_from_charmm):
    bond_pairs = _bond_data_from_charmm.copy()
    for resn, bonds_p in bond_pairs.items():
        bonds_p = bonds_p.split()
        bonds_p = np.array(list(filter(lambda c: c!="BOND", bonds_p)))
        bonds_p = bonds_p.reshape([-1, 2])
        bonds_p = np.array([[atom1, atom2] for atom1, atom2 in bonds_p if not("H" in atom1 or "H" in atom2)])
        bond_pairs[resn] = bonds_p
    angle_pairs = {}
    for resn, bonds_p in bond_pairs.items():
        angles_p = []
        for a1, a2 in bonds_p:
            for b1, b2 in bonds_p:
                if a1==b1 and a2==b2: continue
                if   a1==b1: angle_pair = [a2, b2]
                elif a1==b2: angle_pair = [a2, b1]
                elif a2==b1: angle_pair = [a1, b2]
                elif a2==b2: angle_pair = [a1, b1]
                else: continue
                angles_p.append(angle_pair)
        angle_pairs[resn] = np.array(angles_p)
    return bond_pairs, angle_pairs

def isin_pair(df_k_all, exist_bonds, is_sort=False):
    enm_pairs = np.array(df_k_all[['Atomname_i','Atomname_j']])
    isin_pairs = []
    for p1, p2 in enm_pairs:
        p1_is_exist_bonds = np.isin(exist_bonds, p1)
        p2_is_exist_bonds = np.isin(exist_bonds, p2)
        if is_sort:
            is_exist_bonds = np.append(p1_is_exist_bonds[:, [0]], p2_is_exist_bonds[:, [1]], axis=1)
        else:
            is_exist_bonds = np.logical_or(p1_is_exist_bonds, p2_is_exist_bonds)
        is_exist_bonds = np.any(np.all(is_exist_bonds, axis=1))
        isin_pairs.append(is_exist_bonds)
    isin_pairs = np.array(isin_pairs)
    return isin_pairs

def check_atomname_ij_is_bond(df, sequence, bond_pairs):
    masks_bond = []
    for resid, resn_one_letter in enumerate(sequence):
        if resid==0: continue
        mask_bond1 = df['Resid_i']==resid
        mask_bond2 = df['Resid_j']==resid
        mask_bond3 = isin_pair(df, bond_pairs[f'D{resn_one_letter}'])
        mask_bondn = np.logical_and.reduce([mask_bond1, mask_bond2, mask_bond3])
        masks_bond.append(mask_bondn)
    mask_bond = np.logical_or.reduce(masks_bond)
    ## check lose pair since BOND O3' +3P
    any_resn, next_res_sym = "DG", '+'
    mask_next_res = np.vectorize(lambda s: next_res_sym in s)(bond_pairs[any_resn]).any(axis=1)
    next_res_pairs = bond_pairs[any_resn][mask_next_res]
    next_res_pairs = np.char.replace(next_res_pairs, next_res_sym, '')
    mask_bond_next_res1 = isin_pair(df, next_res_pairs)
    mask_bond_next_res2 = df['Resid_j'] - df['Resid_i'] == 1
    mask_bond_next_res = mask_bond_next_res1 & mask_bond_next_res2
    return mask_bond | mask_bond_next_res

def check_atomname_ij_is_angle(df, sequence, angle_pairs):
    masks_angle = []
    for resid, resn_one_letter in enumerate(sequence):
        if resid==0: continue
        mask_angle1 = df['Resid_i']==resid
        mask_angle2 = df['Resid_j']==resid
        mask_angle3 = isin_pair(df, angle_pairs[f'D{resn_one_letter}'])
        mask_anglen = np.logical_and.reduce([mask_angle1, mask_angle2, mask_angle3])
        masks_angle.append(mask_anglen)
    mask_angle = np.logical_or.reduce(masks_angle)
    ## check lose pair since BOND O3' +P
    next_res_pairs = [["C3'","P"],["O3'","O5'"],["O3'","O1P"],["O3'","O2P"]]
    mask_angle_next_res1 = isin_pair(df, next_res_pairs, is_sort=True)
    mask_angle_next_res2 = df['Resid_j'] - df['Resid_i'] == 1
    mask_angle_next_res = mask_angle_next_res1 & mask_angle_next_res2
    return mask_angle | mask_angle_next_res

def check_ssdna_atomname_ij_is_bond_or_angle(df, sequence, bond_pairs, angle_pairs):
    return np.logical_or(
        check_atomname_ij_is_bond(df, sequence, bond_pairs),
        check_atomname_ij_is_angle(df, sequence, angle_pairs)
    )

def check_dsdna_atomname_ij_is_bond_or_angle(df, sequence, bond_pairs, angle_pairs):
    mask_strand1 = (df['Strand_i']=='STRAND1') & (df['Strand_j']=='STRAND1')
    mask_strand2 = (df['Strand_i']=='STRAND2') & (df['Strand_j']=='STRAND2')
    mask_strand1_is_ba = check_ssdna_atomname_ij_is_bond_or_angle(df[mask_strand1], sequence["STRAND1"], bond_pairs, angle_pairs)
    mask_strand2_is_ba = check_ssdna_atomname_ij_is_bond_or_angle(df[mask_strand2], sequence["STRAND2"], bond_pairs, angle_pairs)
    mask_is_ba = (mask_strand1 & mask_strand1_is_ba) | (mask_strand2 & mask_strand2_is_ba)
    return mask_is_ba

def GetG4DNACategoriesMask(df, df_tetrad_geometry, sequence): ## 參考論文表B.4
    bond_pairs, angle_pairs = GetBondAnglePair()
    mask_0l = np.logical_not(np.logical_or(df['Strand_i']=='STRAND3', df['Strand_j']=='STRAND3'))
    mask_1l = np.logical_xor(df['Strand_i']=='STRAND3', df['Strand_j']=='STRAND3')
    mask_2l = np.logical_and(df['Strand_i']=='STRAND3', df['Strand_j']=='STRAND3')
    mask_0m = check_atomname_ij_is_moiety(df, "M", 0, 'g4dna')
    mask_1m = check_atomname_ij_is_moiety(df, "M", 1, 'g4dna')
    mask_2m = check_atomname_ij_is_moiety(df, "M", 2, 'g4dna')
    mask_0b = check_atomname_ij_is_moiety(df, "B", 0, 'g4dna')
    mask_1b = check_atomname_ij_is_moiety(df, "B", 1, 'g4dna')
    mask_2b = check_atomname_ij_is_moiety(df, "B", 2, 'g4dna')
    # mask_0p = check_atomname_ij_is_moiety(df, "P", 0, 'g4dna')
    # mask_1p = check_atomname_ij_is_moiety(df, "P", 1, 'g4dna')
    # mask_2p = check_atomname_ij_is_moiety(df, "P", 2, 'g4dna')
    mask_0g = check_resid_ij_is_nucleobase(df, "G", 0, sequence['STRAND1'])
    mask_1g = check_resid_ij_is_nucleobase(df, "G", 1, sequence['STRAND1'])
    mask_2g = check_resid_ij_is_nucleobase(df, "G", 2, sequence['STRAND1'])
    mask_same_resid = df['Resid_i']==df['Resid_j']
    mask_diff_resid = np.logical_not(mask_same_resid)
    mask_same_layer = check_resid_ij_same_layer(df, df_tetrad_geometry)
    mask_diff_layer = np.logical_not(mask_same_layer)
    mask_same_strand = check_resid_ij_same_strand(df, df_tetrad_geometry)
    mask_diff_strand = np.logical_not(mask_same_strand)
    mask_is_ba = check_ssdna_atomname_ij_is_bond_or_angle(df, sequence['STRAND1'], bond_pairs, angle_pairs)
    mask_no_ba = np.logical_not(mask_is_ba)
    m2mask = { 
    "dd"   : mask_2l,
    "dg"   : mask_1l,
    "mm"   : mask_0l & mask_2m,
    "bm"   : mask_0l & mask_1m,
    "bbb"  : mask_0l & mask_0m & mask_2b & mask_same_resid & mask_is_ba,
    "bb0"  : mask_0l & mask_0m & mask_2b & mask_same_resid & mask_no_ba,
    "hb"   : mask_0l & mask_0m & mask_2b & mask_diff_resid & mask_2g & mask_same_layer,
    "ss_st": mask_0l & mask_0m & mask_2b & mask_diff_resid & mask_2g & mask_diff_layer & mask_same_strand,
    "cs_st": mask_0l & mask_0m & mask_2b & mask_diff_resid & mask_2g & mask_diff_layer & mask_diff_strand,
    "cl_st": mask_0l & mask_0m & mask_2b & mask_diff_resid & mask_1g,
    "ll_st": mask_0l & mask_0m & mask_2b & mask_diff_resid & mask_0g,
    "brb"  : mask_0l & mask_0m & mask_1b & mask_is_ba,
    "rb"   : mask_0l & mask_0m & mask_1b & mask_no_ba,
    "bpr0" : mask_0l & mask_0m & mask_0b & mask_same_resid & mask_is_ba,
    "pr0"  : mask_0l & mask_0m & mask_0b & mask_same_resid & mask_no_ba,
    "bpr1" : mask_0l & mask_0m & mask_0b & mask_diff_resid & mask_is_ba,
    "pr1"  : mask_0l & mask_0m & mask_0b & mask_diff_resid & mask_no_ba,
    }
    return m2mask

def GetDsDNACategoriesMask(df, sequence, n_bp): ## 參考論文表B.4
    bond_pairs, angle_pairs = GetBondAnglePair()
    mask_0b = check_atomname_ij_is_moiety(df, "B", 0, 'dsdna')
    mask_1b = check_atomname_ij_is_moiety(df, "B", 1, 'dsdna')
    mask_2b = check_atomname_ij_is_moiety(df, "B", 2, 'dsdna')
    # mask_0p = check_atomname_ij_is_moiety(df, "P", 0, 'dsdna')
    # mask_1p = check_atomname_ij_is_moiety(df, "P", 1, 'dsdna')
    # mask_2p = check_atomname_ij_is_moiety(df, "P", 2, 'dsdna')
    mask_same_strand = df['Strand_i']==df['Strand_j']
    mask_diff_strand = np.logical_not(mask_same_strand)
    mask_same_resid = df['Resid_i']==df['Resid_j']
    mask_diff_resid = np.logical_not(mask_same_resid)
    mask_same_resid_1 = (df['Resid_i']-df['Resid_j']==1) | (df['Resid_i']-df['Resid_j']==-1)
    mask_diff_resid_1 = np.logical_not(mask_same_resid_1)
    mask_same_comp_resid = df['Resid_i']+df['Resid_j']==n_bp+1
    mask_diff_comp_resid = np.logical_not(mask_same_comp_resid)
    mask_same_comp_resid_1 = (df['Resid_i']+df['Resid_j']==n_bp) | (df['Resid_i']+df['Resid_j']==n_bp+2)
    mask_diff_comp_resid_1 = np.logical_not(mask_same_comp_resid_1)
    mask_is_ba = check_dsdna_atomname_ij_is_bond_or_angle(df, sequence, bond_pairs, angle_pairs)
    mask_no_ba = np.logical_not(mask_is_ba)
    m2mask = {
    "bbb"  : mask_2b & mask_same_strand & mask_same_resid & mask_is_ba,
    "bb0"  : mask_2b & mask_same_strand & mask_same_resid & mask_no_ba,
    "ss_st": mask_2b & mask_same_strand & mask_diff_resid & mask_same_resid_1,
    "other1":mask_2b & mask_same_strand & mask_diff_resid & mask_diff_resid_1,
    "hb"   : mask_2b & mask_diff_strand & mask_same_comp_resid,
    "cs_st": mask_2b & mask_diff_strand & mask_diff_comp_resid & mask_same_comp_resid_1,
    "other2":mask_2b & mask_diff_strand & mask_diff_comp_resid & mask_diff_comp_resid_1,
    "brb"  : mask_1b & mask_is_ba,
    "rb"   : mask_1b & mask_no_ba,
    "bpr0" : mask_0b & mask_same_resid & mask_is_ba,
    "pr0"  : mask_0b & mask_same_resid & mask_no_ba,
    "bpr1" : mask_0b & mask_diff_resid & mask_is_ba,
    "pr1"  : mask_0b & mask_diff_resid & mask_no_ba,
    }
    return m2mask



def GetCategories(category2mask):
    categories = np.zeros_like(list(category2mask.values())[0], dtype='str')
    for category, mask in category2mask.items():
        cat = np.where(mask, category, '')
        categories = np.char.add(categories, cat)
    return categories
