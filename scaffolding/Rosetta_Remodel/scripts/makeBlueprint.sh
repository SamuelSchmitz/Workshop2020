#!/bin/bash

cat 1gwq_0001.pdb | tail -255 | head -251 | awk '{print $1}' > temp.blueprint
sed -i 's/:NtermProteinFull//g' temp.blueprint 
sed -i 's/:CtermProteinFull//g' temp.blueprint 
sed -i 's/ALA/A/g' temp.blueprint
sed -i 's/CYS/C/g' temp.blueprint
sed -i 's/ASP/D/g' temp.blueprint
sed -i 's/GLU/E/g' temp.blueprint
sed -i 's/PHE/F/g' temp.blueprint
sed -i 's/GLY/G/g' temp.blueprint
sed -i 's/HIS/H/g' temp.blueprint
sed -i 's/ILE/I/g' temp.blueprint
sed -i 's/LYS/K/g' temp.blueprint
sed -i 's/LEU/L/g' temp.blueprint
sed -i 's/MET/M/g' temp.blueprint
sed -i 's/ASN/N/g' temp.blueprint
sed -i 's/PRO/P/g' temp.blueprint
sed -i 's/GLN/Q/g' temp.blueprint
sed -i 's/ARG/R/g' temp.blueprint
sed -i 's/SER/S/g' temp.blueprint
sed -i 's/THR/T/g' temp.blueprint
sed -i 's/VAL/V/g' temp.blueprint
sed -i 's/TRP/W/g' temp.blueprint
sed -i 's/TYR/Y/g' temp.blueprint

sed -i 's/_/\ /g' temp.blueprint
cat temp.blueprint | awk '{print $2 " " $1 " ."}' > 1gwqB.blueprint
rm temp.blueprint
