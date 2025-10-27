from .multiWindows import TimeAgent, SmallTimeAgent
from .makeDataframe import SpSmallTimeAgent
from . import mda
from . import pairtype

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from IPython.display import display, HTML

class DataFrameAgent(SpSmallTimeAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        SpSmallTimeAgent.__init__(self, beg_ns, end_ns, framesperns)
        self.df_all_k = self.read_k_b0_pairtype_df()
        
    def SeleDataFrame(self, df, atomid_i, atomid_j, order=False):
        mask1 = (df['Atomid_i']==atomid_i) & (df['Atomid_j']==atomid_j)
        mask2 = (df['Atomid_i']==atomid_j) & (df['Atomid_j']==atomid_j)
        if order:
            return df[mask1]
        else:
            return df[mask1 | mask2]

    def GetIdentifyPair(self, mode, vlv_threshold=0.1, pe_threshold=1e11, filter=True):
        eigenvector_mul_eigenvector_t = self.eigenvector_mul_eigenvector_t_with_every_mode[mode]
        argsort_i, argsort_j = np.unravel_index(np.argsort(-np.abs(eigenvector_mul_eigenvector_t), axis=None), eigenvector_mul_eigenvector_t.shape)
        df_vv = []
        for i, (ai, aj) in enumerate(zip(argsort_i, argsort_j)):
            vv = eigenvector_mul_eigenvector_t[ai, aj]
            vlv = vv * self.get_eigenvalue_by_id(mode)
            if np.abs(vv)<vlv_threshold: break
            atomid_i = self.df_node.iloc[ai].atomid
            atomid_j = self.df_node.iloc[aj].atomid
            sele = self.SeleDataFrame(self.df_m,  atomid_i, atomid_j, order=filter).squeeze()
            if ai==aj:
                if filter: continue
                pair_id = '-'
                k = '-'
                precentage_error = -1
            elif sele.k.shape!=(): ## k has 0 or >1 values
                if filter: continue
                pair_id = 'x'
                k = 'x'
                precentage_error = -1
            else:
                pair_id = sele.name
                k = sele.k
                precentage_error = np.abs(vlv-k)/(k+1e-9)
            if filter and ((precentage_error < 0) or (pe_threshold < precentage_error)): continue
            df_vv.append([pair_id, sele.Strand_i, sele.Resid_i, sele.Atomname_i, sele.Strand_j, sele.Resid_j, sele.Atomname_j, vv, vlv, k, precentage_error])
        df_vv = pd.DataFrame(df_vv, columns=['PairID', 'Strand_i', 'Resid_i','Atomname_i', 'Strand_j', 'Resid_j','Atomname_j', 'vv', 'vlv', 'k', 'pe'])
        df_vv = df_vv.set_index('PairID')
        return df_vv
        
    def get_hotspot_atom(self, sequence, mode, strandid):
        eigenvector = self.get_eigenvector_by_id(mode)
        idx1 = np.argsort(np.abs(eigenvector))[-1]
        idx2 = np.argsort(np.abs(eigenvector))[-2]
        df_idx1 = self.df_node.iloc[idx1]
        df_idx2 = self.df_node.iloc[idx2]
        hotspot_df_idx = df_idx1 if df_idx1.resid<df_idx2.resid else df_idx2
        hotspot_atom   = f'{hotspot_df_idx.atomname}({hotspot_df_idx.resname[0]}{hotspot_df_idx.resid})'
        if strandid is None and hotspot_df_idx.strand=='STRAND2':
            hotspot_resid = self.n_bp - hotspot_df_idx.resid + 1
        else:
            hotspot_resid = hotspot_df_idx.resid
        return hotspot_atom, hotspot_resid
        
    def transfer_info_from_dict2str(self, mode_info):
        atom_pairs = [value for key, value in mode_info.items() if type(key)==int]
        str_atom_pairs = '\n'.join(atom_pairs)
        # out = f"mode={mode_info['mode']}\nλ={mode_info['lambda']}\n{str_atom_pairs}"
        out = r'$\rm{\lambda}^{\rm{'+self.m_abbr+'}}_{'+str(mode_info['mode'])+r'}$='+str(mode_info['lambda'])\
            +r', ${\{\rm{a}_i\rm{a}_j\}}^{\rm{'+self.m_abbr+r'}}'+'_{'+str(mode_info['mode'])+'}'+r'$:' +'\n'+str_atom_pairs
        return out
        
    def get_modes_info(self, strandid=None, n_mode=80, vlv_threshold=0.1, pe_threshold=0.5, is_return=True):
        dict_modes_info = {}
        dict_resids_info = {}
        for sele_mode in tqdm(range(1, n_mode+1)):
            if not strandid is None:
                modes_in_strandid = np.array(range(1, len(self.eigv_strandid)+1))[self.eigv_strandid==strandid]
                mode_count_by_strandid = modes_in_strandid[sele_mode-1]
            else:
                mode_count_by_strandid = sele_mode
            identify_pair = self.GetIdentifyPair(mode_count_by_strandid, vlv_threshold=vlv_threshold, pe_threshold=pe_threshold, filter=True)
            eigenvalue = self.get_eigenvalue_by_id(mode_count_by_strandid)
            hotspot_atom, _ = self.get_hotspot_atom(self.strandid2sequence, mode_count_by_strandid, strandid) ## TODO: delete hotspot residue       
            if not strandid is None:
                dict_mode_info = {'mode':sele_mode, 'strand':strandid, 'origin_mode': mode_count_by_strandid, 'lambda': np.around(eigenvalue, 1), 'hotspot_atom': hotspot_atom }
            else:
                dict_mode_info = {'mode':sele_mode, 'lambda': np.around(eigenvalue, 1), 'hotspot_atom': hotspot_atom }
            pair_id = 0
            resid_i, resid_j = None, None
            hotspot_resid_i_resid_j = None
            for i,(_, df_iter) in  enumerate(identify_pair.iterrows()):
                resid_i, resid_j = df_iter['Resid_i'], df_iter['Resid_j']
                atomname_i, atomname_j = df_iter['Atomname_i'], df_iter['Atomname_j']
                strand_i, strand_j = df_iter['Strand_i'], df_iter['Strand_j']
                resname_i, resname_j = self.strandid2sequence[strand_i][resid_i], self.strandid2sequence[strand_j][resid_j]                    
                if self.type_na!='dsdna' and self.m_abbr in ['st', 'hb']:
                    if self.get_is_switch_ij(resid_i, resid_j):
                        resid_i, resid_j       = resid_j, resid_i
                        atomname_i, atomname_j = atomname_j, atomname_i
                        strand_i, strand_j     = strand_j, strand_i
                        resname_i, resname_j   = resname_j, resname_i
                if i==0:
                    hotspot_resid_i_resid_j = self.get_hotspot_resid_i_resid_j(resid_i, resid_j)
                atom_pair = f'{atomname_i}{atomname_j}({resname_i}{resid_i}{resname_j}{resid_j})'
                dict_mode_info[pair_id] = atom_pair
                pair_id += 1
            str_mode_info = self.transfer_info_from_dict2str(dict_mode_info)
            dict_modes_info[sele_mode] = dict_mode_info
            if not hotspot_resid_i_resid_j is None:
                dict_resids_info.setdefault(hotspot_resid_i_resid_j, {})
                dict_resids_info[hotspot_resid_i_resid_j][len(dict_resids_info[hotspot_resid_i_resid_j])] = str_mode_info
        df_info_by_modes = pd.DataFrame(dict_modes_info).fillna('').T
        df_info_by_resid = pd.DataFrame(dict_resids_info).fillna('').T
        df_info_by_resid.index.name = 'Residue Pair'
        index = sorted(df_info_by_resid.index, key=sort_key_residue_pair)
        df_info_by_resid = df_info_by_resid.loc[index]
        if is_return:
            return df_info_by_modes, df_info_by_resid
        else:
            self.df_info_by_modes = df_info_by_modes
            self.df_info_by_resid = df_info_by_resid

    def get_layer_and_gstrand_by_resid(self, resid):
        df_gstrand_level_resid = self.df_tetrad_geometry.stack().reset_index()
        row_col = df_gstrand_level_resid[df_gstrand_level_resid[0] == resid][['level_0', 'level_1']].iloc[0]
        return row_col['level_0'], row_col['level_1']
    
    def display_df_info_by_resid(self, save=False):
        display( HTML( self.df_info_by_resid.to_html().replace("\\n","<br>")))
        if save:
            df_info_by_resid_file = f'{self.result_folder}/{self.system}_{self.__class__.__name__}.csv'
            print(f'Write {df_info_by_resid_file}')
            self.df_info_by_resid.to_csv(df_info_by_resid_file)


def sort_key_residue_pair(residue_pair):
    import re
    if not residue_pair:  ## Empty string should always be last
        return (float('inf'), float('inf'), float('inf'))
    elif not ',' in residue_pair: ## sequence position: T5, T5T6
        residue_pair_split = re.findall(r"([A-Za-z]+)(\d+)", residue_pair) ## ex: 'G20G21' -> ['G','20','G','21']
        resid = int(residue_pair_split[0][1])
        return (resid, float('inf'), float('inf'))
    else: ## space position  Q3,top 
        priority_layer   = {'top': 0, 'mid': 1, 'bot': 2}
        priority_gstrand = {'Q3': 0, 'Q2': 1, 'Q1': 2, 'Q4': 3}
        priority_st_type = {'ss_st': 0, 'cs_st': 1, '': 2}

        gstrand, layer = residue_pair.split(',')  ## Split into 'QxQy' and 'top'
        gstrand = gstrand.split('\n') ## Split 'ss_st'
        if len(gstrand)>1: 
            st_stype, gstrand = gstrand
        else:
            st_stype, gstrand = '', gstrand[0]
        gstrand = gstrand[:2]                     ## only Qx
        return (
            priority_st_type[st_stype], ## Sort by st_type: ss_st, cs_st
            priority_layer[layer],      ## Sort by layer: top, mid, bot
            priority_gstrand[gstrand],  ## Sort by the priority of Qx 
        )    

class GraphAgent(DataFrameAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        DataFrameAgent.__init__(self, beg_ns, end_ns, framesperns)
        self.u, self.u_enm = self.load_mda()
        self.df_atom = self.build_df_atom()
        self.df_m = self.get_df_m()
        
    def load_mda(self):
        u = mda.Universe(self.inpcrd_file)
        u_enm = mda.Universe(self.enmcrd_file)
        return u, u_enm

    def build_df_atom(self):
        df_atoms = pd.DataFrame({
        'atomid'   : self.u.atoms.ids,
        'enm'      : self.u_enm.atoms.names,
        'strand'   : self.u.atoms.segids,
        'resid'    : self.u.atoms.resids,
        'resname'  : self.u.atoms.resnames,
        'atomname' : self.u.atoms.names,
        })
        return df_atoms
            
    def spectral_decomposition(self):
        self.build_df_node()
        self.initialize_three_mat()
        self.set_adjacency_by_df(self.df_m)
        self.make_adjacency_symmetry()
        self.build_degree_from_adjacency()
        self.build_laplacian_by_adjacency_degree()
        self.eigen_decompose()
    
    def initialize_three_mat(self):
        self.adjacency_mat = np.zeros((self.n_node, self.n_node))
        self.degree_mat = np.zeros((self.n_node, self.n_node))
        self.laplacian_mat = np.zeros((self.n_node, self.n_node))
        self.b0_mat = np.zeros((self.n_node, self.n_node))
        print('Initialize adjacency, degree and Laplacian matrices... Done.')

    def set_adjacency_by_df(self, df_sele):
        atomid2idx = dict(zip(self.df_node['atomid'], self.df_node.index))
        for _, (atomid_i, atomid_j, k) in self.df_m[['Atomid_i', 'Atomid_j', 'k']].iterrows():
            idx_i = atomid2idx[atomid_i]
            idx_j = atomid2idx[atomid_j]
            self.adjacency_mat[idx_i, idx_j] = k
        
    def make_adjacency_symmetry(self):
        i_lower = np.tril_indices(self.n_node, -1)
        self.adjacency_mat[i_lower] = self.adjacency_mat.transpose()[i_lower]  # make the matrix symmetric
    
    def build_degree_from_adjacency(self):
        for idx in range(self.n_node):
            self.degree_mat[idx, idx] = self.adjacency_mat[idx, :].sum()

    def build_laplacian_by_adjacency_degree(self):
        self.laplacian_mat = self.degree_mat + self.adjacency_mat
        print("Finish the setup for Laplaican matrix.")

    def eigen_decompose(self):
        w, v = np.linalg.eig(self.laplacian_mat)
        idx = w.argsort()[::-1] # sort from big to small
        self.w = w[idx].real
        self.v = v[:, idx].real

    def get_eigenvalue_by_id(self, sele_id):
        return self.w[sele_id-1]

    def get_eigenvector_by_id(self, sele_id):
        return self.v[:,sele_id-1]
        
    def GetEigenvectorEigenvectorT(self):
        eigenvector_mul_eigenvector_t_with_every_mode = {}
        for mode in range(1, self.n_node+1):
            eigenvector = np.array([self.get_eigenvector_by_id(mode)])
            eigenvector_mul_eigenvector_t_with_every_mode[mode] = eigenvector*eigenvector.T
        self.eigenvector_mul_eigenvector_t_with_every_mode = eigenvector_mul_eigenvector_t_with_every_mode

    def get_eigenvector_mode_dict(self, mode=1, _strandid=None):
        ## create empty eigenvector_mode_dict
        atomname_unique = self.df_node.atomname.unique()
        eigenvector_mode_dict = {atomname:{} for atomname in atomname_unique}
        ## generate eigenvector_mode_dict
        eigenvector_mode = self.get_eigenvector_by_id(mode)
        for i, d in self.df_node.iterrows():
            if _strandid == d.strand or _strandid == None:
                eigenvector_mode_dict[d.atomname][d.resid] = eigenvector_mode[i]
        ## dict2list
        for atomname, mode_dict in eigenvector_mode_dict.items():
            eigenvector_mode_dict[atomname] = [mode_dict.get(key, 0.0) for key in range(1, max(mode_dict.keys())+1)]
        return eigenvector_mode_dict

    def plot_identify_atom(self, mode=1, _strandid=None, sign=False):
        self.eigenvector_mode_dict = self.get_eigenvector_mode_dict(mode=mode, _strandid=_strandid)
        num_column = max([len(atomnames) for atomnames in self.atomnamess])
        fig, axes = plt.subplots(len(self.atomnamess), num_column, figsize = (11, 2.2*len(self.atomnamess)))
        for i, atomnames in enumerate(self.atomnamess):
            for j in range(num_column):
                axes[i,j].set_xticks([])
                axes[i,j].get_yaxis().set_visible(False)
            for j, atomname in enumerate(atomnames):
                eigenvector_atomname = self.eigenvector_mode_dict.setdefault(atomname, [])
                if sign:
                    axes[i,j].text(0.6, 0.78, atomname, fontsize=12)
                    axes[i,j].bar(list(range(1, len(eigenvector_atomname)+1)), eigenvector_atomname, color='orange')
                    axes[i,j].set_ylim(-1, 1)
                    axes[i,j].axhline(0.1, linestyle='--')
                    axes[i,j].axhline(0.0, linestyle='--', c='red')
                    axes[i,j].axhline(-0.1, linestyle='--')
                else:
                    axes[i,j].text(0.6, 0.88, atomname, fontsize=12)
                    axes[i,j].bar(list(range(1, len(eigenvector_atomname)+1)), np.abs(eigenvector_atomname), color='orange')
                    axes[i,j].set_ylim(0, 1)
                    axes[i,j].axhline(0.1, linestyle='--')
                if self.m_abbr is 'bm' and atomname is self.get_metal_name_charmm_format():
                    axes[i,j].set_xticks([1, 2])
                else:
                    axes[i,j].set_xticks([3, 9, 15, 21])
                if j==0:
                    axes[i,j].set_ylabel(r'$|v^{st}_{'+str(mode)+'}|$')
                    axes[i,j].get_yaxis().set_visible(True)

        fig.supxlabel('Residue ID', y=0.04)
        plt.subplots_adjust(wspace=0, hspace=0.2)
        plt.show()


class BM(GraphAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        self.m = 'base-metal'
        self.m_abbr = 'bm'
        GraphAgent.__init__(self, beg_ns, end_ns, framesperns)
        self.atomnamess = [
        ## base
        ['C2', 'C4', 'C5', 'C6', 'C7', 'C8'],
        ['N1', 'N2', 'N3', 'N6', 'N7', 'N9'],
        ['O2', 'O4', 'O6'],
        ## metal
        [self.get_metal_name_charmm_format()],
        ]
    def build_df_node(self):  
        atom_type = np.array(list(map(lambda atomname: pairtype.d_atomcgtype.get(atomname), self.df_atom.atomname)))
        mask1 = self.df_atom.strand=='STRAND1'
        mask2 = self.df_atom.strand=='STRAND2'
        mask = mask1 | mask2
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")
    def get_df_m(self):
        self.categories = ['bm']
        mask_1 = np.isin(self.df_all_k['Category'], self.categories)
        df_m = self.df_all_k[mask_1]
        return df_m
    def get_hotspot_resid_i_resid_j(self, resid_i, resid_j):
        metal_name = self.get_metal_name(is_return=True)
        hotspot_resid = f'{metal_name}{resid_j}'
        return hotspot_resid
    def get_metal_name_charmm_format(self):
        mask_mm = self.df_all_k['Category']=='mm'
        metal_name = self.df_all_k[mask_mm]['Atomname_j'].values[0]
        return metal_name
    def get_metal_name(self, is_return=False):
        metal_name_charmm_format = self.get_metal_name_charmm_format()
        if metal_name_charmm_format=='POT':
            metal_name = 'K'
        elif metal_name_charmm_format=='SOD':
            metal_name = 'Na'
        else:
            raise RuntimeError('Could not recongize metal ions', metal_name_charmm_format)
        if is_return:
            return metal_name
        else:
            self.metal_name = metal_name


class HB(GraphAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        self.m = 'hydrogen-bond'
        self.m_abbr = 'hb'
        self.atomnamess = [
        ## base
        ['C2', 'C4', 'C5', 'C6', 'C7', 'C8'],
        ['N1', 'N2', 'N3', 'N6', 'N7', 'N9'],
        ['O2', 'O4', 'O6'],
        ]
        GraphAgent.__init__(self, beg_ns, end_ns, framesperns)
    def build_df_node(self):  
        atom_type = np.array(list(map(lambda atomname: pairtype.d_atomcgtype.get(atomname), self.df_atom.atomname)))
        mask1 = self.df_atom.strand!='STRAND3'
        mask2 = atom_type=='B'
        if self.type_na=='dsdna':
            mask3 = self.df_atom.resid==self.df_atom.resid
        else:
            mask3 = np.isin(self.df_atom.resid, self.df_tetrad_geometry.values.flatten())
        mask = mask1 & mask2 & mask3
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")
    def get_df_m(self):
        self.categories = ['hb']
        mask1 = np.isin(self.df_all_k['Category'], self.categories)
        mask   = mask1
        df_m = self.df_all_k[mask]
        return df_m
    def get_hotspot_resid_i_resid_j(self, resid_i, resid_j):
        if self.type_na == 'dsdna':
            resname_i, resname_j = self.strandid2sequence['STRAND1'][resid_i], self.strandid2sequence['STRAND1'][resid_j]
            hotspot_resid = f'{resname_i}{resid_i}'
        else:
            layer_i, gstrand_i = self.get_layer_and_gstrand_by_resid(resid_i)
            layer_j, gstrand_j = self.get_layer_and_gstrand_by_resid(resid_j)
            hotspot_resid = f'{gstrand_i}{gstrand_j},{layer_i}'
        return hotspot_resid
    def get_is_switch_ij(self, resid_i, resid_j):
            layer_i, gstrand_i = self.get_layer_and_gstrand_by_resid(resid_i)
            layer_j, gstrand_j = self.get_layer_and_gstrand_by_resid(resid_j)
            layer2tetrad_polarity = {
                'top': ['Q3Q2','Q2Q1','Q1Q4','Q4Q3'] if self.GBAc[2]=='anti' else ['Q3Q4','Q4Q1','Q1Q2','Q2Q3'],
                'mid': ['Q3Q2','Q2Q1','Q1Q4','Q4Q3'] if self.GBAc[3]=='anti' else ['Q3Q4','Q4Q1','Q1Q2','Q2Q3'],
                'bot': ['Q3Q2','Q2Q1','Q1Q4','Q4Q3'] if self.GBAc[4]=='anti' else ['Q3Q4','Q4Q1','Q1Q2','Q2Q3'],
            }
            return not gstrand_i+gstrand_j in layer2tetrad_polarity[layer_i]
    

class ST(GraphAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        self.m = 'stack'
        self.m_abbr = 'st'
        self.atomnamess = [
        ## base
        ['C2', 'C4', 'C5', 'C6', 'C7', 'C8'],
        ['N1', 'N2', 'N3', 'N6', 'N7', 'N9'],
        ['O2', 'O4', 'O6'],
        ]
        GraphAgent.__init__(self, beg_ns, end_ns, framesperns)
    def build_df_node(self):  
        atom_type = np.array(list(map(lambda atomname: pairtype.d_atomcgtype.get(atomname), self.df_atom.atomname)))
        mask1 = self.df_atom.strand!='STRAND3'
        mask2 = atom_type=='B'
        mask3 = np.isin(self.df_atom.resid, self.df_tetrad_geometry.values.flatten())
        mask = mask1 & mask2 & mask3
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")
    
class ST_dsdna(ST):
    def get_df_m(self):
        self.categories = ['ss_st']
        mask1 = (self.df_all_k['Strand_i']=='STRAND1') & (self.df_all_k['Strand_j']=='STRAND1')
        mask2  = np.isin(self.df_all_k['Category'], self.categories)
        mask   = mask1 & mask2
        df_m = self.df_all_k[mask]
        return df_m
    def build_df_node(self):  
        atom_type = np.array(list(map(lambda atomname: pairtype.d_atomcgtype.get(atomname), self.df_atom.atomname)))
        mask1 = self.df_atom.strand=='STRAND1'
        mask2 = atom_type=='B'
        mask = mask1 & mask2
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")
    def get_hotspot_resid_i_resid_j(self, resid_i, resid_j):
        resname_i, resname_j = self.strandid2sequence['STRAND1'][resid_i], self.strandid2sequence['STRAND1'][resid_j]
        hotspot_resid = f'{resname_i}{resid_i}{resname_j}{resid_j}'
        return hotspot_resid

class ST_g4dna(ST):
    def get_df_m(self):
        self.categories = ['ss_st', 'cs_st']
        mask_1 = np.isin(self.df_all_k['Category'], self.categories)
        mask_2 = self.df_all_k['Strand_i']=='STRAND1'
        mask_3 = self.df_all_k['Strand_j']=='STRAND1'
        mask   = mask_1 & mask_2 & mask_3
        df_m = self.df_all_k[mask]
        return df_m
    def get_hotspot_resid_i_resid_j(self, resid_i, resid_j):
        if np.abs(resid_i - resid_j) == 1:
            layer_i, gstrand_i = self.get_layer_and_gstrand_by_resid(resid_i)
            layer_j, gstrand_j = self.get_layer_and_gstrand_by_resid(resid_j)
            layer = layer_i if layer_i!='mid' else layer_j
            hotspot_resid = f'ss_st\n{gstrand_i},{layer}'
        else:
            layer_i, gstrand_i = self.get_layer_and_gstrand_by_resid(resid_i)
            layer_j, gstrand_j = self.get_layer_and_gstrand_by_resid(resid_j)
            layer = layer_i if layer_i!='mid' else layer_j
            hotspot_resid = f'cs_st\n{gstrand_i}{gstrand_j},{layer}'
        return hotspot_resid
    def get_is_switch_ij(self, resid_i, resid_j, select_layer=None):
        if select_layer is None:
            select_layer = self.select_layer
        if select_layer=='top':
            return not resid_i in self.df_tetrad_geometry.T['top'].values
        elif select_layer=='bot':
            return not resid_i in self.df_tetrad_geometry.T['mid'].values


class ST_TOP(ST_g4dna):
    def get_df_m(self):
        self.categories = ['ss_st', 'cs_st']
        self.select_layer = 'top'
        mask_1 = np.isin(self.df_all_k['Category'], self.categories)
        mask_i = np.isin(self.df_all_k['Resid_i'], self.df_tetrad_geometry.loc[self.select_layer])
        mask_j = np.isin(self.df_all_k['Resid_j'], self.df_tetrad_geometry.loc[self.select_layer])
        mask   = mask_1 & (mask_i | mask_j)
        df_m = self.df_all_k[mask]
        return df_m
            
class ST_BOT(ST_g4dna):
    def get_df_m(self):
        self.categories = ['ss_st', 'cs_st']
        self.select_layer = 'bot'
        mask_1 = np.isin(self.df_all_k['Category'], self.categories)
        mask_i = np.isin(self.df_all_k['Resid_i'], self.df_tetrad_geometry.loc[self.select_layer])
        mask_j = np.isin(self.df_all_k['Resid_j'], self.df_tetrad_geometry.loc[self.select_layer])
        mask   = mask_1 & (mask_i | mask_j)
        df_m = self.df_all_k[mask]
        return df_m


class RB(GraphAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        self.m = 'ribose-base'
        self.m_abbr = 'rb'
        self.atomnamess = [
        ## base
        ['C2', 'C4', 'C5', 'C6', 'C7', 'C8'],
        ['N1', 'N2', 'N3', 'N6', 'N7', 'N9'],
        ['O2', 'O4', 'O6'],
        ## ribose
        ["C1'", "C2'", "C3'", "C4'", "C5'"],
        ["O3'", "O4'", "O5'"],
        ## phosphate
        ['P', 'O1P', 'O2P']
        ]
        GraphAgent.__init__(self, beg_ns, end_ns, framesperns)
    def build_df_node(self):  
        atom_type = np.array(list(map(lambda atomname: pairtype.d_atomcgtype.get(atomname), self.df_atom.atomname)))
        mask1 = self.df_atom.strand=='STRAND1'
        mask2 = (atom_type=='P') | (atom_type=='R') | (atom_type=='B')
        mask = mask1 & mask2
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")
    def get_df_m(self):
        self.categories = ['rb']
        mask1 = (self.df_all_k['Strand_i']=='STRAND1') & (self.df_all_k['Strand_j']=='STRAND1')
        mask2 = np.isin(self.df_all_k['Category'], self.categories)
        mask   = mask1 & mask2
        df_m = self.df_all_k[mask]
        return df_m
    def get_hotspot_resid_i_resid_j(self, resid_i, resid_j):
        resname_i, resname_j = self.strandid2sequence['STRAND1'][resid_i], self.strandid2sequence['STRAND1'][resid_j]
        hotspot_resid = f'{resname_i}{resid_i}'
        return hotspot_resid


class PR(GraphAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        self.m = 'phosphate-ribose'
        self.m_abbr = 'pr'
        self.atomnamess = [
        ## ribose
        ["C1'", "C2'", "C3'", "C4'", "C5'"],
        ["O3'", "O4'", "O5'"],
        ## phosphate
        ['P', 'O1P', 'O2P']
        ]
        GraphAgent.__init__(self, beg_ns, end_ns, framesperns)
    def build_df_node(self):  
        atom_type = np.array(list(map(lambda atomname: pairtype.d_atomcgtype.get(atomname), self.df_atom.atomname)))
        mask1 = self.df_atom.strand=='STRAND1'
        mask2 = (atom_type=='P') | (atom_type=='R')
        mask  = mask1 & mask2
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        # self.df_node.loc[:, 'idx'] = self.df_node.index
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")

class PR0(PR):
    def __init__(self, beg_ns, end_ns, framesperns):
        super().__init__(beg_ns, end_ns, framesperns)
        self.m_abbr = 'pr0'
    def get_df_m(self):
        self.categories = ['pr0']
        mask1 = (self.df_all_k['Strand_i']=='STRAND1') & (self.df_all_k['Strand_j']=='STRAND1')
        mask2 = np.isin(self.df_all_k['Category'], self.categories)
        mask  = mask1 & mask2
        df_m = self.df_all_k[mask]
        return df_m
    def get_hotspot_resid_i_resid_j(self, resid_i, resid_j):
        resname_i, resname_j = self.strandid2sequence['STRAND1'][resid_i], self.strandid2sequence['STRAND1'][resid_j]
        hotspot_resid = f'{resname_i}{resid_i}'
        return hotspot_resid
    
class PR1(PR):
    def __init__(self, beg_ns, end_ns, framesperns):
        super().__init__(beg_ns, end_ns, framesperns)
        self.m_abbr = 'pr1'
    def get_df_m(self):
        self.categories = ['pr1']
        mask1 = (self.df_all_k['Strand_i']=='STRAND1') & (self.df_all_k['Strand_j']=='STRAND1')
        mask2 = np.isin(self.df_all_k['Category'], self.categories)
        mask  = mask1 & mask2
        df_m = self.df_all_k[mask]
        return df_m
    def get_hotspot_resid_i_resid_j(self, resid_i, resid_j):
        resname_i, resname_j = self.strandid2sequence['STRAND1'][resid_i], self.strandid2sequence['STRAND1'][resid_j]
        hotspot_resid = f'{resname_i}{resid_i}{resname_j}{resid_j}'
        return hotspot_resid

    
class DG(GraphAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        self.m = 'drug_g4dna'
        self.m_abbr = 'dg'
        GraphAgent.__init__(self, beg_ns, end_ns, framesperns)
    def build_df_node(self):  
        mask = self.df_atom==self.df_atom
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")
    def get_df_m(self):
        self.categories = ['dg']
        mask_1 = np.isin(self.df_all_k['Category'], self.categories)
        df_m = self.df_all_k[mask_1]
        return df_m
    
class CLST(GraphAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        self.m = 'core-loop-stack-or-hbond'
        self.m_abbr = 'cl_st'
        GraphAgent.__init__(self, beg_ns, end_ns, framesperns)
    def build_df_node(self):  
        mask = self.df_atom==self.df_atom
        self.df_node = self.df_atom[mask].reset_index(drop=True)
        self.n_node = len(self.df_node)
        print(f"Thare are {self.n_node} nodes.")
    def get_df_m(self):
        self.categories = ['cl_st']
        mask_1 = np.isin(self.df_all_k['Category'], self.categories)
        df_m = self.df_all_k[mask_1]
        return df_m