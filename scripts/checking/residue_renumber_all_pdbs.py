#! /usr/bin/env python

#This script renumbers all pdb files that are in PDB path
# residue_renumber_all_pdbs.py <PDB path> <GROMACS programs path>
# Example: python residue_renumber_all_pdbs.py /home/faccioli/Execute/1VII/ /home/faccioli/Programs/gmx-4.6.5/no_mpi/bin/

import sys
import os
from subprocess import Popen, PIPE

def run_renumber_by_editconf(pdbfile, gmx_path):
	program = os.path.join(gmx_path, "editconf")
	process = Popen([program, '-f', pdbfile, '-o', pdbfile, '-resnr', '1'], stdout=PIPE, stderr=PIPE)
	stdout, stderr = process.communicate()	 	

def get_files_pdb(mypath):
	only_tar_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".pdb"):
				f_path = os.path.join(root,file)
				only_tar_file.append(f_path)			
	return only_tar_file

def main():
	#Avoid GROMACS backup files
	os.environ["GMX_MAXBACKUP"]="-1"

	#path of pdb files
	path_pdb = sys.argv[1]
	#GROMACS program path
	gmx_path = sys.argv[2] 

	for pdbfile in get_files_pdb(path_pdb):
		run_renumber_by_editconf(pdbfile, gmx_path)

main()