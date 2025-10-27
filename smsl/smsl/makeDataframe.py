from .enmBuilder import ENMAgent
from .trajTransfer import MDAgent
from .multiWindows import TimeAgent, SmallTimeAgent
from .category import GetDsDNACategoriesMask, GetG4DNACategoriesMask, GetCategories

import pandas as pd
import numpy as np

class SpringAgent(ENMAgent, MDAgent):
    def __init__(self):
        super().__init__()
    def load_TimeAgent(self):
        self.t_agent = TimeAgent(SpSmallTimeAgent, create_mean=True)
        self.t_agent.mean.initial_folder(verbose=1)
        self.t_agent.copy_mean_crd()
    ## writing csv file
    def write_all_k_b0_pairtype_df(self):
        atoms = self.u_noh.atoms
        enm_names = self.u_enm.atoms.names
        ## get atom info
        enm_name2atomid = {enm_name:i for i, enm_name in enumerate(enm_names, start=1)}
        atomid2strand = {i:resid for i, resid in enumerate(atoms.chainIDs, start=1)}
        atomid2resid = {i:resid for i, resid in enumerate(atoms.resids, start=1)}
        atomid2atomname = {i: atomname for i, atomname in enumerate(atoms.names, start=1)}
        for time_label, st_agent in self.t_agent.items():
            st_agent.df_all_k = st_agent.generate_k_b0_pairtype_df(enm_name2atomid=enm_name2atomid, atomid2strand=atomid2strand, atomid2resid=atomid2resid, atomid2atomname=atomid2atomname)
        st_agent = list(self.t_agent.values())[0]
        ## get category
        if self.type_na == 'g4dna':
            get_category = GetG4DNACategoriesMask 
            catgory2mask = get_category(st_agent.df_all_k, self.df_tetrad_geometry, self.strandid2sequence) 
        elif self.type_na == 'dsdna':
            get_category = GetDsDNACategoriesMask
            catgory2mask = get_category(st_agent.df_all_k, self.strandid2sequence, self.n_bp) 
        else:
            print('could not determine the type_na')
        for time_label, st_agent in self.t_agent.items():
            st_agent.df_all_k['Category'] = GetCategories(catgory2mask)
            st_agent.write_k_b0_pairtype_df()
        ## get mean k value
        self.t_agent.mean.df_all_k = self.calc_k_mean()
        self.t_agent.mean.write_k_b0_pairtype_df()
    def calc_k_mean(self):
        ks_per_windows = []
        b0_per_windows = []
        for time_label, st_agent in self.t_agent.items():
            ks_per_windows.append(st_agent.df_all_k['k'].values)
            b0_per_windows.append(st_agent.df_all_k['b0'].values)
        df_all_k = st_agent.df_all_k.copy()
        df_all_k['k'] = np.array(ks_per_windows).mean(axis=0)
        df_all_k['b0'] = np.array(b0_per_windows).mean(axis=0)
        return df_all_k
    ## reading csv file
    def read_all_k_b0_pairtype_df(self):
        for time_label, st_agent in self.t_agent.items():
            st_agent.df_all_k = st_agent.read_k_b0_pairtype_df()
        self.t_agent.mean.df_all_k = self.t_agent.mean.read_k_b0_pairtype_df()


class SpSmallTimeAgent(SmallTimeAgent):
    def generate_k_b0_pairtype_df(self, **kwargs):
        l_enm_i = np.genfromtxt(self.end_prm_file, skip_header=4, skip_footer=1, dtype=str, usecols=[0]).T
        l_enm_j = np.genfromtxt(self.end_prm_file, skip_header=4, skip_footer=1, dtype=str, usecols=[1]).T
        l_k     = np.genfromtxt(self.end_prm_file, skip_header=4, skip_footer=1, dtype=float, usecols=[2]).T
        l_b0    = np.genfromtxt(self.end_prm_file, skip_header=4, skip_footer=1, dtype=float, usecols=[3]).T
        l_pair_id  = np.arange(1, len(l_k)+1)
        l_aotmid_i = np.vectorize(kwargs['enm_name2atomid'].get)(l_enm_i)
        l_aotmid_j = np.vectorize(kwargs['enm_name2atomid'].get)(l_enm_j)
        l_strand_i = np.vectorize(kwargs['atomid2strand'].get)(l_aotmid_i)
        l_strand_j = np.vectorize(kwargs['atomid2strand'].get)(l_aotmid_j)
        l_resid_i  = np.vectorize(kwargs['atomid2resid'].get)(l_aotmid_i)
        l_resid_j  = np.vectorize(kwargs['atomid2resid'].get)(l_aotmid_j)
        l_atomname_i = np.vectorize(kwargs['atomid2atomname'].get)(l_aotmid_i)
        l_atomname_j = np.vectorize(kwargs['atomid2atomname'].get)(l_aotmid_j)
        df_all_k = pd.DataFrame({
            'PairID'    : l_pair_id,
            'Enm_i'     : l_enm_i,
            'Strand_i'  : l_strand_i,
            'Resid_i'   : l_resid_i,
            'Atomname_i': l_atomname_i,
            'Atomid_i'  : l_aotmid_i,
            'Enm_j'     : l_enm_j,
            'Strand_j'  : l_strand_j,
            'Resid_j'   : l_resid_j,
            'Atomname_j': l_atomname_j,
            'Atomid_j'  : l_aotmid_j,
            'k'         : l_k,
            'b0'        : l_b0,
        })
        return df_all_k
    def write_k_b0_pairtype_df(self):
        self.df_all_k.to_csv(self.df_all_k_file, index=False)
        print('Writing data to:', self.df_all_k_file)
    def read_k_b0_pairtype_df(self):
        print('Reading data from:', self.df_all_k_file)
        return pd.read_csv(self.df_all_k_file, index_col='PairID')