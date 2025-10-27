from . import mda
from .multiWindows import TimeAgent, SmallTimeAgent
from .config import ConfAgent


class ENMAgent(ConfAgent):
    def __init__(self):
        super().__init__()
    def load_TimeAgent(self, is_return=False):
        t_agent = TimeAgent(SmallTimeAgent, create_mean=False)
        t_agent.initialTimeAgent_folder(verbose=0)
        if is_return:
            return t_agent
        else:
            self.t_agent = t_agent
    def load_ENMUniverse(self, t_agent=None, is_return=False):
        if t_agent is None: t_agent = self.t_agent
        all_crd_file = [tg.inpcrd_file for tg in t_agent.values()]
        npt2_crd_file = all_crd_file[0]
        u_enm = mda.Universe(npt2_crd_file, *all_crd_file)
        if is_return:
            return u_enm
        else:
            print(str(u_enm.trajectory).split('more ')[-1][:-1])
            self.u_enm = u_enm
    def modify_ENMUniverse(self, u_enm=None):
        if u_enm is None: u_enm = self.u_enm
        ## modify resname
        u_enm.residues.resnames = 'NA'
        u_enm.residues.resids = 1
        ## modify atomname
        enm_chain, enm_id = '@', 0
        tmp_chainID = ''
        for atom_enm in u_enm.atoms:
            ## generate A1, A2 ...
            if atom_enm.segid != tmp_chainID:
                tmp_chainID = atom_enm.segid
                enm_chain = chr(ord(enm_chain)+1)
                enm_id = 1
            else:
                enm_id += 1
            atom_enm.name = f'{enm_chain}{enm_id}'
            atom_enm.type = atom_enm.name
    def write_all_crd_file(self, t_agent=None, u_enm=None):
        if t_agent is None: t_agent=self.t_agent
        if u_enm is None: u_enm = self.u_enm
        for i, (time_label, tg) in enumerate(t_agent.items()):
            u_enm.atoms.write(tg.enmcrd_file, frames=[i])
    def get_union_bonds(self, t_agent=None, u_enm=None):
        if t_agent is None: t_agent=self.t_agent
        if u_enm is None: u_enm = self.u_enm
        ## Get the union of all d_pairs in every windows
        for _, time_label in zip(u_enm.trajectory, t_agent.keys()):
            bonds_per_frame = []
            for atom_enm in u_enm.atoms:
                atoms_pair = u_enm.select_atoms(f'around {self.cutoff} id {atom_enm.id}')
                bonds_per_atom = [atom_enm+atom_pair for atom_pair in atoms_pair]
                bonds_per_frame.extend(bonds_per_atom)
            print(f'{time_label}ns: {len(bonds_per_frame)//2} pairs')
            u_enm.add_bonds(bonds_per_frame)
        print(f'\nTotal {len(u_enm.atoms)} beads with {len(u_enm.bonds)} pairs')
    def write_all_crd_rtf_str_psf_file(self, t_agent=None, u_enm=None):
        if t_agent is None: t_agent=self.t_agent
        if u_enm is None: u_enm = self.u_enm
        rtfOut = self.MakeEnmRTF()
        strOut = self.MakeEnmStr()
        for i, (time_label, tg) in enumerate(t_agent.items()):
            print(f'{time_label}ns')
            u_enm.atoms.write(tg.enmcrd_file, frames=[i])
            self.PostModifyEnmCrd(i, tg)
            open(tg.enmstr_file, "w").write(self.strOut)
            open(tg.enmrtf_file, "w").write(self.rtfOut)
            u_enm.atoms.convert_to.parmed().write_psf(tg.enmpsf_file)
    def PostModifyEnmCrd(self, i, tg):
        with open(tg.enmcrd_file, 'r+') as file:
            content = file.read()  # Read entire content
            if content.count('*')==2: ## if not modify yet
                modified_content = content.replace(f"{i} FROM", "FROM\n*")  # Modify the content
                file.seek(0)  # Move to the beginning of the file
                file.write(modified_content)  # Write the modified content back to the file
                file.truncate()  # If the new content is shorter than the old, truncate remaining
    def MakeEnmStr(self, u_enm=None, is_return=False, split_count=299):
        if u_enm is None: u_enm = self.u_enm
        strOut, bond_count = f'* {self.system}\n', 0
        for bond in u_enm.bonds:
            if bond_count==0:
                strOut += 'IC EDIT\n'
            atom1 = bond.atoms[0].name
            atom2 = bond.atoms[1].name
            strOut += f'DIST 1 {atom1:5} 1 {atom2:5} 0.0\n'
            bond_count += 1
            if bond_count>=split_count:
                strOut += 'END\n'
                bond_count = 0
        if not strOut.endswith('END\n'):
            strOut += 'END\n'
        if is_return:
            return strOut
        else:
            self.strOut = strOut
    def MakeEnmRTF(self, u_enm=None, is_return=False):
        if u_enm is None: u_enm = self.u_enm
        rtfOut = ''
        ## 1. write title
        rtfOut += f'''\
* {self.system}
*
41  1
'''
        ## 2. write mass
        for a in u_enm.atoms:
            rtfOut += f'MASS{a.id:6} {a.name:<6} {a.mass:8.6f}\n'
        ## 3. write defa_res
        rtfOut += f'''\

DEFAULT FIRS NONE LAST NONE
RESI NA 0.00
GROUP
'''
        ## 4. write atoms
        for a in u_enm.atoms:
            rtfOut += f'ATOM {a.name:<5}{a.name:<5}   0.00\n'
        ## 5. write bonds
        for b in u_enm.bonds:
            rtfOut += f'BOND {b.atoms[0].name:5} {b.atoms[1].name:<5}\n'
        ## 5. write end
        rtfOut += '\n\n\nEND'
        if is_return:
            return rtfOut
        else:
            self.rtfOut = rtfOut