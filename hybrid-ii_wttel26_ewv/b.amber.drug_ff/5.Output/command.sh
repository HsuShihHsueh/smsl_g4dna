pdbid='2jpz'
drug='EWV'
complex=${pdbid}'_'${drug}

cp ../1.GenTopology/${drug}.amb2gmx/*${drug}* ./
cp ../4.Align/${drug}_GMX.align.gro.pdb ./${drug}_GMX.gro.pdb
cp ../4.Align/${pdbid}.pgro.pdb ./

## a. Add `TER` between `DNA (TER) ION` maunally
sed -i '/DNA/{
    n;/ION/i TER
}' ${pdbid}.pgro.pdb

## b.合併G四聯體與藥物
cat ${pdbid}.gro.pdb    > ${complex}.gro.pdb ## ${pdbid}.gro.pdb 由 gmx pdb2gmx -f ${pdbid}.pgro.pdb  -o ${pdbid}.gro.pdb (../c.gromacs.run_md/1.create_system/1.1_create_system.ipynb)而來，更換成force field 的 atomname
cat ${drug}_GMX.gro.pdb >> ${complex}.gro.pdb
# Add 3 `TER` between `DNA (TER) K (TER) Drug (TER) END` maunally

## c.Modify ${drug}_GMX.top
# 1. 將[nbfunc, Compound, ]關鍵字的下一行註解
sed -i '/nbfunc/{N;s/\n/\n;/}' ${drug}_GMX.top
sed -i '/Compound/{N;s/\n/\n;/}' ${drug}_GMX.top
# 2. 將 [atomtypes] 內 name bond_type的文字轉小寫
sed -E '/\[ atomtypes \]/,/^$/ {
    /^[[:space:]]*[A-Z]/ {
        s/^([[:space:]]*)([A-Z0-9*]+)([[:space:]]+)([A-Z0-9*]+)/\1\L\2\3\L\4/
    }
}' ${drug}_GMX.top > temp.top && mv temp.top ${drug}_GMX.top
# 3. 將 [atoms] 內 type的文字轉小寫
sed -E '/\[ atoms \]/,/^$/ {
    /^[[:space:]]*[0-9]/ {
        s/^([[:space:]]*[0-9]+[[:space:]]+)([A-Z0-9*]+)([[:space:]]+)/\1\L\2\3/
    }
}' ${drug}_GMX.top > temp.top && mv temp.top ${drug}_GMX.top



