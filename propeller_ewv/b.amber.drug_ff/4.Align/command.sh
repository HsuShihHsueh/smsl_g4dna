## step1. align 6ccw(hybrid-ii_ewv, mobile) to 1kf1.pgro.pdb(propeller, target) via pymol
## and save the ewv
pymol >>
load 1kf1.pgro.pdb ## from propeller_ewv/a.charmm.dna_setup/1kf1.pgro.pdb
load 6ccw.pdb ## from https://files.rcsb.org/download/6ccw.pdb
align 6ccw and resi 4-6+22-24, 1kf1.pgro and resi 2-4+20-22
set retain_order, 1
save 6ccw_EWV_align_1kf1.pdb, 6ccw and resn EWV
>> EOF

## step2. rename atomname maunally 
## 6ccw_EWV_align_1kf1.pdb -> 6ccw_EWV_align_1kf1_reatomname.pdb

## step3. EWV_GMX.gro(optimize by g16) -> EWV_GMX.gro.pdb
gmx editconf -f ../2.Mdgx/EWV.amb2gmx/EWV_GMX.gro -o EWV_GMX.gro.pdb

## step4. align ewv(optimize by g16, mobile) to ewv(step2, target) via pymol
pymol >>
load EWV_GMX.gro.pdb ## from step3
load 6ccw_EWV_align_1kf1_reatomname.pdb ##　from step2
align EWV_GMX.gro, 6ccw_EWV_align_1kf1_reatomname
set retain_order, 1
save EWV_GMX.align.gro.pdb, EWV_GMX.gro
>> EOF


## Since EWV has align to 6ccw in system:hybrid-ii_ewv, 
## so just copy ../../hybrid-ii_ewv/b.amber.drug_ff/4.Align/EWV_GMX.gro.pdb

