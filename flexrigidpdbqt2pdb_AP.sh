#!/bin/bash

# flexrigidpdbqt2pdb.sh_revAP

# ProSciens, Computing & Molecular Sciences
# http://www.prosciens.com/

# Edelmiro Moman (2011) 

# This script is distributed under the GNU General Public License (GPL v3.0)
# http://www.gnu.org/licenses/gpl.html

## This script merges the rigid and the flexible PDBQT output files produced by AutoDock4 and AudoDock Vina and converts the recombined file to PDB format. It namely:

# 1. Extracts one solution from a combined flexible PDBQT file (or a single-model PDBQT when MODEL is set to 0)
# 2. Converts the extracted model to PDB format by using either pdbqt_to_pdb.py (DEFAULT) or OpenBabel (version >= 2.3)
# 3. Converts the rigid PDBQT file to PDB format by using either pdbqt_to_pdb.py (DEFAULT) or OpenBabel (version >= 2.3)
# 4. Merges both PDB files
# 5. Produces both the apo protein and (when LIGAND=1) the holo protein-ligand complex
# 6. Moves the output files to a given directory (when an OUTDIR is provided)

# The script assumes that both the rigid and the flexible files have a .PDBQT extension
# The script assumes that OpenBabel is installed, it is in the user's path and it is invoked by the command "babel"
# If BABEL is set to 1, MGLTools is not required to run the script. However, the version of OpenBabel should the 2.3 or higher in order to handle PDBQT files
# If BABEL is set to 0, both MGLTools and OpenBabel are required, but any recent version of OpenBabel should work
# For more information on the required software, please, visit:
# http://autodock.scripps.edu/
# http://openbabel.org/wiki/Main_Page
# This far only AutoDock Vina 1.1.1, MGLTools-1.5.4 and OpenBabel 2.3 for GNU/Linux have been tested
# You are welcome to post bug reports as well as improvements (prosciens.com) 

###################################################################################
#                                                                                 #
#                                 OPTIONS                                         #
#                                                                                 # 
#  The user may need/want to configure variables 1-8 before running the script:   #
#                                                                                 #
###################################################################################

# 1. Name of the rigid PDBQT file without extension

RIGIDPDBQT=$1

# 2. Name of the output PDBQT without extension

FLEXPDBQT=$2

# 3. Name of the ligand

LIGNAME=$3

# 4. Number of the model to be extracted from the flexible PDBQT solutions. By default, the first model will be extracted. If the file contains a single model set MODEL=0 (use only if absolutely certain)

MODEL=1

# 5. If you have MGLTools installed, please, enter the full installation path

MGLTOOLS=/home/andrea/Programs/mgltools_x86_64Linux2_1.5.7rc1

# 6. Otherwise, we will use OpenBabel for the conversion (you need OpenBabel 2.3 or higher, it needs to be in your path and to be called "babel") 
# "0" means "no" and "1" means "yes"

BABEL=0

# 7. Do you want to have the ligand in the merged PDB file? 
# "0" means "no" and "1" means "yes"

LIGAND=1

# 8. Do you want to copy the merged file(s) to a specific directory? If yes, please, provide the path, otherwise, just leave blank

OUTDIR=$4

###################################
#                                 #
#       The fun starts here       #
#                                 #
###################################

###  Let's extract the model

if [ "$MODEL" == "0" ]
then
cp $FLEXPDBQT.pdbqt ${FLEXPDBQT}_${MODEL}.pdbqt
else
sed -n '/MODEL '$MODEL'$/{:a;n;/ENDMDL/b;p;ba}' $FLEXPDBQT.pdbqt > ${FLEXPDBQT}_${MODEL}.pdbqt
fi

###  Let's convert both rigid and flexible PDBQTs to PDBs

if [ "$BABEL" == "1" ]
then
obabel -ipdbqt $RIGIDPDBQT.pdbqt -opdb $RIGIDPDBQT.pdb
obabel -ipdbqt ${FLEXPDBQT}_${MODEL}.pdbqt -opdb ${FLEXPDBQT}_${MODEL}.pdb
else
$MGLTOOLS/bin/pythonsh $MGLTOOLS/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f $RIGIDPDBQT.pdbqt
$MGLTOOLS/bin/pythonsh $MGLTOOLS/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f ${FLEXPDBQT}_${MODEL}.pdbqt
fi

cp $RIGIDPDBQT.pdb ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb

###  Let's merge the files

# 1. First we clean up the model PDB and we divide the flex res and the ligand

grep ATOM ${FLEXPDBQT}_${MODEL}.pdb > ${FLEXPDBQT}_${MODEL}.pdb.tmp
grep HETATM ${FLEXPDBQT}_${MODEL}.pdb > ${LIGNAME}_${MODEL}.pdb.tmp

# 2. Next we create a list of residues

cut -c 18-27 ${FLEXPDBQT}_${MODEL}.pdb.tmp > residuelistraw.tmp
cat residuelistraw.tmp | uniq > residuelist.tmp


# 3. Then we split the model file into residues

while read r
do
rns=`echo $r | sed 's/ //g'`
egrep "[ \t]$r[ \t]" ${FLEXPDBQT}_${MODEL}.pdb.tmp > $rns.pdb.tmp
sed -i 's/'$FLEXPDBQT'_'$MODEL'.pdb.tmp://' $rns.pdb.tmp

# 4. And we do the actual merging

sed '/'"HN  ${r}"'/ r '$rns'.pdb.tmp' ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb > ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb.tmp
mv ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb.tmp ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb
done < residuelist.tmp

# 5. Some atom renumbering and PDB cleanup

obabel -ipdb ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb -opdb ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb

# 6. Merging the ligand

if [ "$LIGAND" == "1" ]
then
##sed -i 's/ATOM  /HETATM/g' LIG.pdb.tmp
grep ATOM ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb > ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb.tmp
cat ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb.tmp ${LIGNAME}_${MODEL}.pdb.tmp > ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_holo.pdb

# 7. Plus some more atom renumbering and cleanup

obabel -ipdb ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_holo.pdb -opdb ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_holo.pdb
fi

###  Optionaly, we can copy the merged file(s) to a given directory

if [ -n "$OUTDIR" ]
then
cp ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_apo.pdb $OUTDIR/
fi
if [[ -n "$OUTDIR" && $LIGAND == "1" ]]
then
cp ${RIGIDPDBQT}_${LIGNAME}_${MODEL}_holo.pdb $OUTDIR/
fi

###  Finally, we clean after ourselves

rm *.tmp

