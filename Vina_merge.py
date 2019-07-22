#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:26:39 2019

@author: Andrea Pasquadibisceglie
"""
from io import StringIO
import subprocess
import argparse
import os


###SOME PATHS###
#mgltools directory path
MGLTOOLSDIR = "/home/user/mgltools_x86_64Linux2_1.5.7rc1/"

#path for utilities and pythonsh of mgltools. In linux OS these should be the same
utilities24 = "{0}MGLToolsPckgs/AutoDockTools/Utilities24/".format(MGLTOOLSDIR)
pythonsh = '{0}bin/pythonsh'.format(MGLTOOLSDIR)


#command line parsing#
parser = argparse.ArgumentParser(prog='Vina_merge', usage='%(prog)s [options]')

Input = parser.add_argument_group('Input')
Input.add_argument('-w', '--worDir', action='store', help='working directory (default: vina output directory)', metavar='')
Input.add_argument('-r', '--rigid', action='store', help='rigid PDBQT', required=True, metavar='')
Input.add_argument('-p', '--poses', action='store', help='vina output PDBQT', required=True, metavar='')
Input.add_argument('-n', '--number',  action='store', help='rank of the docking pose to merge (default: 1)', metavar='')

Output = parser.add_argument_group('Output')
Output.add_argument('-o', '--outDir', action='store', help='output directory (default: vina output directory)', metavar='')

args = parser.parse_args()

#REQUIRED VALUES#
rigid=args.rigid
poses=args.poses

#DEFAULT VALUES#
if args.worDir is None:
    worDir=os.path.dirname(args.poses)
else:
    worDir=args.worDir

if args.number is None:
    number="1"
else:
    number=args.number

if args.outDir is None:
    outDir=worDir
else:
    outDir=args.outDir

#SOME FUNCTIONS#
def pdbqt2pdb(pdbqt):
    subprocess.run([pythonsh, utilities24+"pdbqt_to_pdb.py", "-f {0}".format(pdbqt), "-o {0}".format(pdbqt)[:-2]])


#########################################
###(1)EXTRACT THE MODEL YOU WANT MERGE###
with open (poses, "r") as fi:
    lines = fi.readlines()
    for i, l in enumerate(lines):
        if l.split()[0] == "MODEL" and l.split()[1] == number:
            model_n = l.split()[0] + l.split()[1]
            j=i #remember starting position
            #extract pose n
            with open (worDir+model_n+".pdbqt", "w") as model:
                for k in range(j, len(lines)):
                    if lines[k].startswith("ENDMDL") == False:
                        model.write(lines[k])
                    else:
                        model.write(lines[k])
                        break

###(2)CONVERT PDBQT TO PDB###
pdbqt2pdb(model.name) #io.TextIOWrapper
pdbqt2pdb(rigid)
#save path
modelpdb = (model.name)[:-2]
rigidpdb = (rigid)[:-2]

###(3)SPLIT LIGAND AND FLEX RESIDUES###
with open (modelpdb, "r") as fi:
    flex = StringIO()
    lig = StringIO()
    for lines in fi.readlines():
        if lines.startswith("ATOM"):
            flex.write(lines)
        elif lines.startswith("HETATM"):
            lig.write(lines)

###(4)MAKE LIST OF FLEX RESIDUES###
reslist=[]
contflex = flex.getvalue()
listaflex = contflex.split("\n")
for line in listaflex:
    if line != "":
        reslist.append(line[17:26])

###(5)MERGE###
merge=StringIO()
with open (rigidpdb, "r") as fi:
    for liner in fi.readlines():
        merge.write(liner) #line of the rigid
        if liner[17:26] in reslist:
            for linef in listaflex:
                if linef[17:26] == liner[17:26]:
                    merge.write(linef+"\n") #all lines of the flex residues matched
                    reslist.remove(liner[17:26]) #remove the flex from the list
        if liner.startswith("TER"):
            merge.write(lig.getvalue()) #merge the ligand at the end

###(6)RENUMBERING###
atnum=1
remerge=StringIO()
for line in (merge.getvalue()).split("\n"):
    if line.startswith("ATOM") and int(line.split()[1])==(atnum+1):
        atnum = atnum+1
        remerge.write(line+"\n")
    #check if the atom number is correct
    elif line.startswith("ATOM") and int(line.split()[1])!=(atnum+1):
        #replace the substring using 5 position aligned on left
        newline = line.replace("{0:>5}".format(line.split()[1]),
                               "{0:>5}".format(str(atnum+1)))
        remerge.write(newline+"\n")
        atnum = atnum+1
    else:
        remerge.write(line+"\n")

###(7)WRITE MERGED FILE###
with open (outDir+"merge.pdb", "w") as fo:
    cont = remerge.getvalue()
    fo.write(cont)
