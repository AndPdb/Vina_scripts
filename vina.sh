#!/bin/sh

VINAPATH=/home/andreap/Programs/AutodockVina/vina
LIG=./ligands/*.pdbqt
OUT=./output

for ligand in $LIG

do
	lig=`basename "$ligand" .pdbqt`
	$VINAPATH --ligand $ligand --config conf.txt --out $OUT/out_$lig.pdbqt --log $OUT/log_$lig.txt
done
