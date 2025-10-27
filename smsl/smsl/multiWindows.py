import os
import numpy as np
import shutil
from .config import ConfAgent


class SmallTimeAgent(ConfAgent):
    def __init__(self, beg_ns, end_ns, framesperns):
        super().__init__()
        self.beg_ns = beg_ns
        self.end_ns = end_ns
        self.framesperns = framesperns
        self.initial_file_variable()
    def initial_file_variable(self):
        self.basename = os.path.basename(os.getcwd())
        self.scratch_folder = f'{self.scratchroot}/{self.basename}'
        self.time_label = f'{self.beg_ns:04}_{self.end_ns:04}'
        self.time_list  = [self.beg_ns, self.end_ns]
        self.frame_list = [int(self.beg_ns*self.framesperns), int(self.end_ns*self.framesperns)]
        self.time_folder = f'{self.bigtraj_folder}/{self.time_label}'
        self.stime_folder= f'{self.scratch_folder}/{self.time_label}' 
        self.md_folder   = f'{self.time_folder}/md_data'
        self.enm_folder  = f'{self.time_folder}/enm_data'
        self.script_folder  = f'{self.time_folder}/script'
        self.result_folder  = f'{self.time_folder}/result'
        self.sbackup_folder = f'{self.stime_folder}/backup'
        self.inpcrd_file = f'{self.md_folder}/{self.system}.nohydrogen.crd'
        self.inpdcd_file = f'{self.md_folder}/{self.time_label}.nohydrogen.dcd'
        self.enmcrd_file = f'{self.enm_folder}/na_enm.crd'
        self.enmrtf_file = f'{self.enm_folder}/na_enm_{self.cutoff:.2f}.rtf'
        self.enmstr_file = f'{self.enm_folder}/na_enm_{self.cutoff:.2f}.str'
        self.enmpsf_file = f'{self.enm_folder}/na_enm_{self.cutoff:.2f}.psf'
        self.sprm_file   = f'{self.stime_folder}/na_enm_{self.cutoff:.2f}.prm'
        self.savg_file   = f'{self.stime_folder}/avg_{self.cutoff:.2f}.ic'
        self.sfluct_file = f'{self.stime_folder}/fluct_{self.cutoff:.2f}.ic'
        self.snmainp_file= f'{self.stime_folder}/nma.inp'
        self.snmadat_file= f'{self.stime_folder}/nma.dat'
        self.nmainp_file = f'{self.script_folder}/nma.inp'
        self.error_file  = f'{self.enm_folder}/error_{self.cutoff:.2f}.txt'
        self.end_prm_file= f'{self.enm_folder}/na_enm_{self.cutoff:.2f}.prm'
        self.iterate_file= f'{self.result_folder}/k_fluct_iter_{self.cutoff:.2f}.csv'
        self.df_all_k_file = f'{self.result_folder}/pairtypes_k_b0_cutoff_{self.cutoff:.2f}.csv'
    def initial_folder(self, verbose=1, sfolder=False):
        initial_folders = [self.bigtraj_folder, self.time_folder, self.md_folder, 
                           self.enm_folder, self.script_folder, self.result_folder]
        if sfolder:
            initial_sfolders = [self.scratch_folder, self.stime_folder, self.sbackup_folder]
            initial_folders = initial_folders + initial_sfolders
        for folder in initial_folders:
            if not os.path.exists(folder):
                if verbose: print(f'mkdir {folder}')
                os.makedirs(folder)

class TimeAgent(dict, ConfAgent):
    def __init__(self, SeleSmallTimeAgent, create_mean=False):
        ConfAgent.__init__(self)
        self.initial_dict_value(SeleSmallTimeAgent, create_mean)
    def initial_dict_value(self, SeleSmallTimeAgent, create_mean):
        self.framesperns = self.frame_num / self.time_num
        beg_ns_list = np.linspace(0, self.time_num, self.split_num+2, dtype='int')[:-2]
        end_ns_list = beg_ns_list + (beg_ns_list[1]-beg_ns_list[0])*2
        for beg_ns, end_ns in zip(beg_ns_list, end_ns_list):
            self[f'{beg_ns:04}_{end_ns:04}'] = SeleSmallTimeAgent(beg_ns, end_ns, self.framesperns)
        if create_mean:
            self.mean = SeleSmallTimeAgent(0, self.time_num, self.framesperns)
    def copy_mean_crd(self):
        ## copy crd file from first window to mean
        src_file = list(self.values())[0].inpcrd_file
        dst_file = self.mean.inpcrd_file
        if not os.path.exists(dst_file):
            print(f'cp {src_file} {dst_file}')
            shutil.copyfile(src_file, dst_file)
        src_file = list(self.values())[0].enmcrd_file
        dst_file = self.mean.enmcrd_file
        if not os.path.exists(dst_file):
            print(f'cp {src_file} {dst_file}')
            shutil.copyfile(src_file, dst_file)
    def initialTimeAgent_folder(t_agent, verbose=1, sfolder=False):
        for time_label, st_agent in t_agent.items():
            st_agent.initial_folder(verbose, sfolder)