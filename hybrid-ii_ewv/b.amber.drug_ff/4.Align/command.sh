## step1. align 6ccw(hybrid-ii+ewv, mobile) to 2jpz(hybrid-ii, target) via pymol
## and save the ewv
pymol >>
load 2jpz.pgro.pdb ## from hybrid-ii_ewv/a.charmm.dna_setup/2jpz.pgro.pdb
load 6ccw.pdb ## from https://files.rcsb.org/download/6ccw.pdb
align 6ccw, 2jpz.pgro
set retain_order, 1
save 6ccw_EWV_align_2jpz.pdb, 6ccw and resn EWV
>> EOF

## step2. rename atomname maunally 
## 6ccw_EWV_align_2jpz.pdb -> 6ccw_EWV_align_2jpz_reatomname.pdb

## step3. EWV_GMX.gro(optimize by g16) -> EWV_GMX.gro.pdb
gmx editconf -f ../2.Mdgx/EWV.amb2gmx/EWV_GMX.gro -o EWV_GMX.gro.pdb

## step4. align ewv(optimize by g16, mobile) to ewv(step2, target) via pymol
pymol >>
load EWV_GMX.gro.pdb ## from step3
load 6ccw_EWV_align_2jpz_reatomname.pdb ##　from step2
align EWV_GMX.gro, 6ccw_EWV_align_2jpz_reatomname
set retain_order, 1
save EWV_GMX.align.gro.pdb, EWV_GMX.gro
>> EOF
