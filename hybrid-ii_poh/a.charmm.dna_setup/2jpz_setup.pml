## 1. load pdb file from rcsb
load 2JPZ.pdb
## 2. modify dna resid
remove resi 1-2 or resi 25-26
alter all,resi=str(int(resi)-2)
## 3. load 2 potassium ions
load K.cif
create K2, K
set_name K, K1
## 4. (optional) optimize the 3D view
set sphere_scale, 0.6
zoom
## 5. move potassium to the central of tetrad
translate [-3.137, 0.213,-0.343], K1
translate [ 0.902,-0.001,-0.020], K2
## 6. define the resid of potassium
alter model K1, resi=23
alter model K2, resi=24
## 7. save the molecule of dna and ion part separately to pdb
save 2jpz_dna.pdb, polymer.nucleic
save 2jpz_ion.pdb, inorganic

