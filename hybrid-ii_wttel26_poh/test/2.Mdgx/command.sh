mkdir conf1 iter1
tleap 0.GenTop.tleap
mdgx  -i 1.GenConformers.fix.py
bash     2.OptOrca.sh
mdgx  -i 3.ExtractEnergies.mdgx
mdgx  -i 4.FitParm.mdgx
## combine tmpyp4.all.modify.frcmod+iter1/tmpyp4.dat=iter1/tmpyp4.final.frcmod by maunally
tleap -f 5.GenNewTop.tleap
acpype -p iter1/tmpyp4.prmtop -x iter1/tmpyp4.rst7
## change tmpyp4 name from 'mpy' to 'POH'
'''
multiscale    14 * 16
multiphysics  14 * 16
allostery      7 * 20 + 4(node009)
crounus       48 + 48 + 40
total         728
'''

