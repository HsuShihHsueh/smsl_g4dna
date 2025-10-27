#!/bin/bash

## transfer the atomname from rcsb format(base on IUPAC-IUB rules) to charmm format
## usage           : rcsb2charmm.sh <pdbfile> > <pdbfile>
## example         : rcsb2charmm.sh 1kf1.pdb > 1kf1.0.pdb
## charmm format   : charmm/toppar/top_all36_na.rtf
## rcsb format     : https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf
## IUPAC-IUB rules : Standard residue abbreviations conform to the IUPAC-IUB rules in J. Biol. Chem. 241, 527, 2491 (1966).


# base
sed "s/ DA/ADE/g" $1 |\
sed "s/ DT/THY/g"    |\
sed "s/ DC/CYT/g"    |\
sed "s/ DG/GUA/g"    |\
# phosphate
sed "s/OP1/O1P/g"    |\
sed "s/OP2/O2P/g"    |\
sed "s/OP3/O3P/g"    |\
sed "s/OP4/O4P/g"    |\
# methyl- group of thymine
sed "s/ C7/C5M/g"    |\
sed "s/H71/H51/g"    |\
sed "s/H72/H52/g"    |\
sed "s/H73/H53/g"    |\
sed "s/H74/H54/g"    |\
# 5' & 3' hydrogen on O5' or O3' 
#   H5'    
#   |
# --C5'--O5'
#   |    |
#   H5'' HO5'
sed "s/HO5'/ H5T/g"  |\
sed "s/HO3'/ H3T/g"  |\
# metal ion
sed "s/K  /POT/"     |\
sed "s/  K/POT/"     |\
sed "s/K  /POT/"     |\
sed "s/NA /SOD/"     |\
sed "s/ NA/SOD/"     |\
sed "s/NA /SOD/"     
	
 
