## step1. align 2hri(propeller+poh, mobile) to 1kf1(propeller, target) via pymol
## and save the poh
pymol >>
load 1kf1.pgro.pdb ## from propeller_poh/a.charmm.dna_setup/1kf1.pgro.pdb
load 2hri.pdb ## from https://files.rcsb.org/download/2hri.pdb
align 2hri, 1kf1.pgro
save 2hri_POH_align_1kf1.pdb, 2jpz and resn POH
>> EOF

## step2. rename atomname by maunally 
## 2hri_POH_align_1kf1.pdb -> 2hri_POH_align_1kf1_reatomname.pdb

## step3. POH_GMX.gro(optimize by g16) -> POH_GMX.gro.pdb
gmx editconf -f ../2.Mdgx/POH.amb2gmx/POH_GMX.gro -o POH_GMX.gro.pdb

## step4. align poh(optimize by g16, mobile) to poh(step2, target) via pymol
pymol >>
load POH_GMX.gro.pdb ## from step3
load 2hri_POH_align_1kf1_reatomname.pdb ##　from step2
align POH_GMX.gro, 2hri_POH_align_1kf1_reatomname
save POH_GMX.align.gro.pdb, POH_GMX.gro
>> EOF
