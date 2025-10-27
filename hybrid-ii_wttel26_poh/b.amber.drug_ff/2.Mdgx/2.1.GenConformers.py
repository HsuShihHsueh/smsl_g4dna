from pymol import cmd
from numpy import random
import os

drug = "POH"
n_conf = 360
cpus = 8    ## using {cpus} cpus in parallel computing 
rams = 7000 ## using {rams} MB in computing
w_text = f'''\
! PAL{cpus}
! mp2 cc-pvdz TightSCF  NoKeepInts

%scf
  MaxCore {rams}
end
%mp2
  MaxCore {rams}
end
'''

cmd.load(f'../1.GenTopology/{drug}.mol2')
cmd.set('retain_order', 1)
cmd.alter('all', f"resn='{drug}'")
cmd.alter('all', "elem=''")
for i in range(1, n_conf+1):
    cmd.set_dihedral('name C12', 'name C11', 'name C3', 'name C4', i-random.rand())
    cmd.save(f'conf1/Conf{i}.pdb')
    os.system(f'obabel conf1/Conf{i}.pdb -o orcainp -O conf1/Conf{i}.orca')
    
    with open(f'conf1/Conf{i}.orca', 'r') as f:
        r_text = f.read()
        r_text = r_text[r_text.find('*'):]

    with open(f'conf1/Conf{i}.orca', 'w') as f:
        f.write(w_text+r_text)
exit()
