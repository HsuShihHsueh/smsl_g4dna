## 1. load pdb file from rcsb
load 1KF1.pdb
## 2. remove water of crystallization
remove solvent
## 3. remove the third potassium
remove resi 26 and name K
## 4. (optional) optimize the 3D view
set sphere_scale, 0.6
zoom
## 5. redefine the resid of potassium
alter resi 24 and name K, resi=23
alter resi 25 and name K, resi=24
## 6. save the molecule of dna and ion part separately to pdb
save 1kf1_dna.pdb, polymer.nucleic
save 1kf1_ion.pdb, inorganic

