import os
import datetime
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from scipy import constants
from shutil import copyfile
from subprocess import check_call
from concurrent.futures import ThreadPoolExecutor

from . import mda
from .multiWindows import TimeAgent, SmallTimeAgent
from .config import ConfAgent
from .opt import charmm_c41b1


def exec_charmm(f_input, f_output, charmmexec=None, verbose=1):
    if charmmexec is None:
        charmmexec = os.path.abspath(charmm_c41b1)
    if verbose:
        print(f'$ {charmmexec} < {f_input} > {f_output}')
    check_call(charmmexec, stdin=open(f_input, 'r'), stdout=open(f_output, 'w+'), shell=True)


class FmAgent(ConfAgent):
    def __init__(self):
        super().__init__()
    def load_TimeAgent(self, is_return=False):
        t_agent = TimeAgent(FmSmallTimeAgent, create_mean=False)
        t_agent.initialTimeAgent_folder(verbose=0, sfolder=True)
        ## initial all temperature relate variable: RT, alpha
        for time_label, st_agent in t_agent.items():
            st_agent.initial_temp_var()
        if is_return:
            return t_agent
        else:
            self.t_agent = t_agent
    def write_all_nmascipt(self, t_agent=None):
        if t_agent is None: t_agent = self.t_agent
        for time_label, st_agent in t_agent.items():
            st_agent.write_nmascript()
    def run_all_fluctmatch(self, max_workers=48):
        for time_label, st_agent in self.t_agent.items():
            st_agent.loop_iterate = tqdm(range(1, self.iter_num+1), desc=f'Window {time_label}')
        with ThreadPoolExecutor(max_workers=max_workers) as p:
            for time_label, st_agent in self.t_agent.items():
                p.submit(st_agent.fluct_match)            


class FmSmallTimeAgent(SmallTimeAgent):
    def initial_temp_var(self):    
        self.RT = self.T * (constants.k * constants.N_A / (constants.calorie * constants.kilo)) ## RT kcal/mol # https://en.wikipedia.org/wiki/KT_(energy)
        self.alpha = self.RT * self.lr
    def make_ic_table_from_md(self, RT, is_retrun=False, precision_ic=4, verbose=0):
        self.u_enm = mda.Universe(self.enmpsf_file, self.inpdcd_file)
        coord = self.u_enm.trajectory.timeseries().astype(np.float64)
        bonds_avg, bonds_fluct = [], []
        bonds = tqdm(self.u_enm.bonds) if verbose else self.u_enm.bonds
        for bond in bonds:
            atom1 = bond[0].ix
            atom2 = bond[1].ix
            ## np.linalg axis may be speed up
            bonds_dist = np.linalg.norm(coord[atom1,:,:]-coord[atom2,:,:], axis=1)
            bonds_avg.append(bonds_dist.mean())
            bonds_fluct.append(bonds_dist.std())
        df_ic_table = pd.DataFrame({
            'atom_i': self.u_enm.bonds.atom1.names,
            'atom_j': self.u_enm.bonds.atom2.names,
            'b'     : bonds_avg,
            'fluct0': bonds_fluct,
        })
        if is_retrun:
            return df_ic_table
        else:
            self.df_ic_table = df_ic_table
    def guess_k_value(self, bonds_avg, bonds_fluct, RT):
        bonds_k = RT / np.square(bonds_fluct)
        return bonds_k
    def calc_pow_m2(self, arr):
        pow_m2 = np.reciprocal(np.square(arr))
        return pow_m2
    def read_ic_data(self, filename, skip_header=5):
        data_ic = np.genfromtxt(filename, skip_header=skip_header, dtype=np.float64)[:, 17]
        return data_ic
    def write_prm_file(self, iternum, df_ic_table=None):
        if df_ic_table is None: df_ic_table = self.df_ic_table
        prmOut = f'''\
* {self.system} iternum={iternum}
*

BONDS
'''
        for idx, d in df_ic_table.iterrows():
            prmOut += f"{d['atom_i']:<6} {d['atom_j']:<6}{d['k']:9.4f}  {d['b']:8.4f}\n"
        prmOut += "\nEND"
        open(self.sprm_file, "w").write(prmOut)
        
    def write_nmascript(self):
        nmainpOut = f'''\
* Author: Shih-Hsueh Hsu
* 

bomlev -2

set rootfolder {self.bigtraj_folder}
set scratchfolder {self.scratch_folder}
set timelabel {self.time_label} 
set nadir   @rootfolder/@timelabel
set scnadir @scratchfolder/@timelabel
set mddir   @nadir/md_data
set enmdir  @nadir/enm_data
set crdfile @enmdir/na_enm.crd
set rtffile @enmdir/na_enm_{self.cutoff:.2f}.rtf
set icstr   @enmdir/na_enm_{self.cutoff:.2f}.str
set prmfile @scnadir/na_enm_{self.cutoff:.2f}.prm
set vibcrd  @scnadir/na_enm.vib_{self.cutoff:.2f}.crd 
set vibic   @scnadir/na_enm_{self.cutoff:.2f}.vib 
set avgic   @scnadir/avg_{self.cutoff:.2f}.ic 
set fluctic @scnadir/fluct_{self.cutoff:.2f}.ic 


set fileu  10 
set temp   {self.T:.1f}

open unit 11 read form name @rtffile
read rtf card unit 11 
close unit 11 

open read unit @fileu form name @prmfile 
read para unit @fileu card 
close unit @fileu 

read sequ card 
* 
 1
NA
generate DNA setup

open read unit @fileu card name @crdfile 
read coor unit @fileu ignore 
coor copy comp 
close unit @fileu 

scalar mass set 1 sele all end

skip all excl bond 
update inbfrq 0 1hbfrq 0 imgfrq 0 
ener 
mini sd nstep 100 
mini abnr nstep 45000 

coor orie rms mass 
scalar wmain copy mass 

ioformat extended 
open write unit @fileu card name @vibcrd 
write coor unit @fileu card 

stream @icstr 

ic fill 
open write unit @fileu card name @avgic 
write ic   unit @fileu card resid 
* Internal coordinate averages 
* 

close unit @fileu 

calc nmode   ?natom * 3 

set nmodes   @nmode 
set type     temp 
set fluctu   20 
set vibu     22 

open write unit @fluctu card name @fluctic 
open write unit @vibu   card name @vibic 

vibran nmode @nmodes 
    diag fini 
    fluc ic @type @temp tfre 0.0 mode 7 thru @nmodes 
    ic save 
    ic write unit @fluctu resid 
    * Internal coordinate fluctuation 
    * 
    write normal card mode 1 thru @nmodes unit @vibu 
end 

stop
'''
        open(self.snmainp_file, 'w').write(nmainpOut)
        open(self.nmainp_file, 'w').write(nmainpOut)

    def fluct_match(self, loop_iterate=None, charmmexec=None, verbose=0):
        if loop_iterate is None:
            if self.loop_iterate is None:
                loop_iterate = range(1, self.iter_num+1)
            else:
                loop_iterate = self.loop_iterate
        ## Get icavg, icfluct from MD
        self.make_ic_table_from_md(self.RT, verbose=verbose)
        self.df_iter_table = self.df_ic_table[['atom_i','atom_j']].copy() ## iter_table
        self.df_iter_table['b0'] = self.df_ic_table['b']                  ## iter_table
        self.df_iter_table['f0'] = self.df_ic_table['fluct0']             ## iter_table
        
        ## Guess K
        self.df_ic_table['k'] = self.guess_k_value(self.df_ic_table['b'], self.df_ic_table['fluct0'], self.RT)
        self.df_ic_table['target'] = self.calc_pow_m2(self.df_ic_table['fluct0'])
        self.df_iter_table['k0'] = self.df_ic_table['k']                  ## iter_table
        self.df_iter_table['target'] = self.df_ic_table['target']         ## iter_table
    
        # Update and Backup PRM
        self.write_prm_file(iternum=0)
        prm_backup = f'{self.sbackup_folder}/na_enm_{self.cutoff:.2f}.{0}.prm'
        copyfile(self.sprm_file, prm_backup)
    
        ## Run Charmm
        exec_charmm(self.snmainp_file, self.snmadat_file, charmmexec, verbose=verbose)
    
        ## Updata K
        icavg = self.read_ic_data(self.savg_file)
        self.df_ic_table['fluct'] = self.read_ic_data(self.sfluct_file)
        self.df_ic_table['k'] = self.guess_k_value(icavg, self.df_ic_table['fluct'], self.RT)
        self.df_iter_table['b.1'] = icavg                                  ## iter_table
        self.df_iter_table['f.1'] = self.df_ic_table['fluct']              ## iter_table
        self.df_iter_table['k_init'] = self.df_ic_table['k']               ## iter_table
        
        ## New error file
        errorOut = f'''\
    Created at {datetime.datetime.now()}
    n_iter  error   
    '''
        open(self.error_file, 'w').write(errorOut)
    
        for i in loop_iterate:
            ## Run Charmm
            exec_charmm(self.snmainp_file, self.snmadat_file, charmmexec, verbose=verbose)
    
            # Read icavg, icfluct
            self.df_ic_table['fluct'] = self.read_ic_data(self.sfluct_file)
            optimized = self.calc_pow_m2(self.df_ic_table['fluct'])
    
            # Update k
            self.df_ic_table['k'] -= self.alpha * (optimized - self.df_ic_table['target'])
            self.df_ic_table.loc[self.df_ic_table['k']<0., 'k'] = 0. # Make all negative values zero
    
            # Update iter
            self.df_iter_table[f'f{i}'] = self.df_ic_table['fluct']           ## iter_table
            self.df_iter_table[f'k{i}'] = self.df_ic_table['k']               ## iter_table
            
            # Calculate RMSD error and Write
            error = self.df_ic_table['fluct0'] - self.df_ic_table['fluct']
            error = np.sqrt(np.square(error).mean())
            open(self.error_file, 'a').write(f'{i:<5} {error:8.4f}\n')
    
            # Update and Backup PRM
            self.write_prm_file(iternum=i)
            prm_backup = f'{self.sbackup_folder}/na_enm_{self.cutoff:.2f}.{i}.prm'
            copyfile(self.sprm_file, prm_backup)
        ## Copy the Final PRM
        copyfile(self.sprm_file, self.end_prm_file)
        self.df_iter_table.to_csv(self.iterate_file) ## iter_table