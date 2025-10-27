## write from 2 package:
## https://github.com/yizaochen/fluctmatch
## https://github.com/yizaochen/enmspring
import os
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from MDAnalysis.topology.guessers import guess_types, guess_masses
from . import mda
from .multiWindows import TimeAgent, SmallTimeAgent
from .config import ConfAgent


class MDAgent(ConfAgent):
    def __init__(self):
        super().__init__()
    def load_TimeAgent(self, is_return=False):
        t_agent = TimeAgent(SmallTimeAgent, create_mean=False)
        t_agent.initialTimeAgent_folder(verbose=0)
        if is_return:
            return t_agent
        else:
            self.t_agent = t_agent
    def load_MDUniverse(self, pdbFile=None, xtcFile=None, is_return=False):
        if pdbFile is None:
            pdbFile = os.path.realpath(glob.glob('*.pdb')[0])
        if xtcFile is None:
            xtcFile = os.path.realpath(glob.glob('*.xtc')[0])
        print(f'Reading PDB file: {pdbFile}')
        print(f'Reading XTC file: {xtcFile}')
        u_traj = mda.Universe(pdbFile, xtcFile)
        if is_return:
            return u_traj
        else:
            self.u_traj = u_traj
    def modify_MDUniverse(self, seg_select, u_traj=None):
        if u_traj is None: u_traj = self.u_traj
        print('Split Segment ...')
        self.split_segment(u_traj, seg_select)
        print('Guess Drug Atomtype ...')
        self.guess_atomtype(u_traj)
        print('Transfer Atomtype from RCSB to Charmm Format ...')
        self.pdb_rscb2charmm_format(u_traj)
    def pdb_rscb2charmm_format(self, u):
        ## modify resname from rcsb pdb format to charmm format
        u.select_atoms('nucleic and resname DA').residues.resnames = 'ADE'
        u.select_atoms('nucleic and resname DT').residues.resnames = 'THY'
        u.select_atoms('nucleic and resname DC').residues.resnames = 'CYT'
        u.select_atoms('nucleic and resname DG').residues.resnames = 'GUA'
        ## modify atomname from rcsb pdb format to charmm format
        ### phosphate
        u.select_atoms('nucleic and name OP1').atoms.names = 'O1P'
        u.select_atoms('nucleic and name OP2').atoms.names = 'O2P'
        u.select_atoms('nucleic and name OP3').atoms.names = 'O3P'
        u.select_atoms('nucleic and name OP4').atoms.names = 'O4P'
        ### methyl- group of thymine
        u.select_atoms('nucleic and resname *T* and name C7 ').atoms.names = 'C5M'
        u.select_atoms('nucleic and resname *T* and name H71').atoms.names = 'H51'
        u.select_atoms('nucleic and resname *T* and name H72').atoms.names = 'H52'
        u.select_atoms('nucleic and resname *T* and name H73').atoms.names = 'H53'
        u.select_atoms('nucleic and resname *T* and name H74').atoms.names = 'H54'
        ### 5' & 3' hydrogen on O5' or O3' 
        ###   H5'    
        ###   |
        ### --C5'--O5'
        ###   |    |
        ###   H5'' HO5'
        u.select_atoms("nucleic and name HO5'").atoms.names = 'H5T'
        u.select_atoms("nucleic and name HO3'").atoms.names = 'H3T'
        ### metal ion
        u.select_atoms('name K ').residues.resnames = 'POT'
        u.select_atoms('name NA').residues.resnames = 'SOD'
        u.select_atoms('name K ').atoms.names = 'POT'
        u.select_atoms('name NA').atoms.names = 'SOD'
    def split_segment(self, u, seg_select):
        for seg, selection in seg_select.items():
            ag = u.select_atoms(selection)
            n_residue = len(ag.residues.resids)
            ag.atoms.chainIDs = seg
            ag.residues.resids = range(1, n_residue+1)   
            print(f'  select {seg} for {len(ag):>4} atoms')
    def guess_atomtype(self, u, verbose=1):
        ag = u.atoms[u.atoms.masses==0]
        print(f'  Guess {len(ag)} atoms')
        ag.atoms.elements = guess_types(ag.atoms.names)
        ag.atoms.types    = guess_types(ag.atoms.names)
        ag.atoms.masses   = guess_masses(ag.atoms.types)
        if verbose:
            print('  Elements after guess:\n ', ag.atoms.elements)
    def get_nohydrogen_ag(self, u_traj=None, is_return=False):
        if u_traj is None: u_traj = self.u_traj
        u_noh = u_traj.select_atoms('not element H')
        if is_return:
            return u_noh
        else:
            print(f'selct {len(u_noh)} heavy atoms')
            self.u_noh = u_noh
    def write_crd_and_dcd_file(self, tg):
        beg_ns, end_ns = tg.time_list
        beg_frame, end_frame = tg.frame_list
        output = f'{beg_ns:4d} → {end_ns:4d} ns'
        self.u_noh.write(tg.inpcrd_file, frames=[beg_frame])
        self.u_noh.write(tg.inpdcd_file, frames=list(range(beg_frame, end_frame)))
        print(output)
    def write_all_crd_and_dcd_file(self, u_noh=None, enable_multithreading=True):
        if u_noh is None: u_noh = self.u_noh
        print(f'''\
Whole Trajectory: 
frame_num : {self.frame_num} frames
time_num  : {self.time_num} ns
Split_num : {self.split_num}

===Time_list===
''')
        if enable_multithreading:
            with ProcessPoolExecutor() as executor:
                futures = [executor.submit(self.write_crd_and_dcd_file, tg) for tg in self.t_agent.values()]
                # Wait for all tasks to complete
                for future in futures:
                    future.result()  # This will re-raise any exceptions encountered
        else:
            for tg in self.t_agent.values():
                self.write_crd_and_dcd_file(tg)