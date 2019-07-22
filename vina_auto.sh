#!/bin/sh

### Script for virtual screening. It shuld be run in the same folder of the receptor and conf.txt. ligands shuld be in a subfolder named ligands, the config file shuld be named conf.txt, and the output folder named output.

VINAPATH=/home/andreap/Programs/AutodockVina/vina
LIG=./ligands/*.pdbqt
OUT=./output

for ligand in $LIG

do
	lig=`basename "$ligand" .pdbqt`
	$VINAPATH --ligand $ligand --config conf.txt --out $OUT/out_$lig.pdbqt --log $OUT/log_$lig.txt
done
