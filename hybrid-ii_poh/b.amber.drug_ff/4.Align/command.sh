## step1. align 1kf1.npt2_rms.pdb(system:propeller_poh, mobile) to 2jpz.pgro.pdb(hybrid-ii, target) via pymol
## and save the ewv
pymol >>
load 2jpz.pgro.pdb ## from hybrid-ii_poh/a.charmm.dna_setup/2jpz.pgro.pdb
load 1kf1.npt2_rms.pdb ## from propeller_poh/c.gromacs.run_md/1.create_system/1kf1.npt2_rms.pdb
align 1kf1.npt2_rms and resi 2-4+20-22, 2jpz.pgro and resi 2-4+20-22
set retain_order, 1
save 1kf1_POH_align_2jpz.pdb, 1kf1.npt2_rms and resn POH
>> EOF

## step2. rename atomname maunally (there is no need, just copy)
## 1kf1_POH_align_2jpz.pdb -> 1kf1_POH_align_2jpz_reatomname.pdb

## step3. POH_GMX.gro(optimize by g16) -> POH_GMX.gro.pdb
gmx editconf -f ../2.Mdgx/POH.amb2gmx/POH_GMX.gro -o POH_GMX.gro.pdb

## step4. align poh(optimize by g16, mobile) to poh(step2, target) via pymol
pymol >>
load POH_GMX.gro.pdb ## from step3
load 1kf1_POH_align_2jpz_reatomname.pdb ##　from step2
align POH_GMX.gro, 1kf1_POH_align_2jpz_reatomname
set retain_order, 1
save POH_GMX.align.gro.pdb, POH_GMX.gro
>> EOF

