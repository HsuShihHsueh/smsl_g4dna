#!/bin/bash

## transfer the atomname from charmm format to rcsb format(base on IUPAC-IUB rules)
## usage           : charmm2rcsb.sh <pdbfile> > <pdbfile>
## example         : charmm2rcsb.sh 1kf1.aa.pdb > 1kf1.pgro.pdb
## charmm format   : charmm/toppar/top_all36_na.rtf
## rcsb format     : https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf
## IUPAC-IUB rules : Standard residue abbreviations conform to the IUPAC-IUB rules in J. Biol. Chem. 241, 527, 2491 (1966).


# base
sed "s/ADE/ DA/g" $1 |\
sed "s/THY/ DT/g"    |\
sed "s/CYT/ DC/g"    |\
sed "s/GUA/ DG/g"    |\
# phosphate
sed "s/O1P/OP1/g"    |\
sed "s/O2P/OP2/g"    |\
sed "s/O3P/OP3/g"    |\
sed "s/O4P/OP4/g"    |\
# methyl- group of thymine
sed "s/C5M/ C7/g"    |\
sed "s/H51/H71/g"    |\
sed "s/H52/H72/g"    |\
sed "s/H53/H73/g"    |\
sed "s/H54/H74/g"    |\
# 5' & 3' hydrogen on O5' or O3' 
#   H5'    
#   |
# --C5'--O5'
#   |    |
#   H5'' HO5'
sed "s/ H5T/HO5'/g"  |\
sed "s/ H3T/HO3'/g"  |\
# metal ion
sed "s/POT/K  /"     |\
sed "s/POT/  K/"     |\
sed "s/POT/K  /"     |\
sed "s/SOD/NA /"     |\
sed "s/SOD/ NA/"     |\
sed "s/SOD/NA /"     
	
 
