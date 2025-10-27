
from pymol import cmd
from numpy.random import rand
from os import system

n_conf = 360
w_text = '''\
! PAL2
! mp2 cc-pvdz TightSCF  NoKeepInts

%scf
  MaxCore 512
end
%mp2
  MaxCore 512
end
'''

cmd.load('tmpyp4.mol2')
cmd.set('retain_order', 1)
cmd.alter('all', "resn='mpy'")
cmd.alter('all', "elem=''")
for i in range(1, n_conf+1):
    cmd.set_dihedral('name C12', 'name C11', 'name C3', 'name C4', i-rand())
    cmd.save(f'conf1/Conf{i}.pdb')
    os.system(f'obabel conf1/Conf{i}.pdb -o orcainp -O conf1/Conf{i}.orca')
    
    with open(f'conf1/Conf{i}.orca', 'r') as f:
        r_text = f.read()
        r_text = r_text[r_text.find('*'):]

    with open(f'conf1/Conf{i}.orca', 'w') as f:
        f.write(w_text+r_text)

exit()


