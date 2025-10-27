import os

script_path = os.path.realpath(__file__)
script_folder = os.path.dirname(script_path)

toppar               = os.path.abspath(f'{script_folder}/../opt/c43b1/toppar')
charmm_c43b1         = os.path.abspath(f'{script_folder}/../opt/c43b1/bin/charmm')
charmm_c41b1         = os.path.abspath(f'{script_folder}/../opt/c41b1_yz/bin/charmm')
rcsb2charmmformat_sh = os.path.abspath(f'{script_folder}/../opt/atomname_transfer/rcsb2charmmformat.sh')
charmm2rcsbformat_sh = os.path.abspath(f'{script_folder}/../opt/atomname_transfer/charmm2rcsbformat.sh')



