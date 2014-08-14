#! /usr/bin/env python

#This script check all pdb files using pdb2gmx
# check_structures_gromacs.py <PDB path> <GROMACS programs path>
# Example: python check_structures_gromacs.py /home/faccioli/Execute/1VII/ /home/faccioli/Programs/gmx-4.6.5/no_mpi/bin/

import sys
import os
from subprocess import Popen, PIPE
import shutil

#Log file of structures NOT accepted by pdb2gmx
log_file_structures_not_accepted_by_pdb2gmx='structures_not_accepted_by_pdb2gmx.log'
def get_files_pdb(mypath):
	only_tar_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".pdb"):
				f_path = os.path.join(root,file)
				only_tar_file.append(f_path)			
	return only_tar_file

def get_all_directories(mypath):
	only_directory = []
	for root, dirs, files in os.walk(mypath):
		for directory in dirs:
			only_directory.append(directory)
	return only_directory

def delete_check_files():		
	if os.path.isfile('check.top'):
		os.remove('check.top')	
	if os.path.isfile('posre.itp'):
		os.remove('posre.itp')

def rename_atoms_structue_by_pdb2gmx(pdbfile, gmx_path, forcefield='amber99sb-ildn'):
	delete_check_files()
	program = os.path.join(gmx_path, "pdb2gmx")
	process = Popen([program, '-f', pdbfile, '-o', pdbfile, '-p', 'check.top', '-water', 'none', '-ff', forcefield], stdout=PIPE, stderr=PIPE)
	stdout, stderr = process.communicate()	 	
	delete_check_files()	

def main():
	#Avoid GROMACS backup files
	os.environ["GMX_MAXBACKUP"]="-1"

	#path of pdb files
	path_pdb = sys.argv[1]
	#GROMACS program path
	gmx_path = sys.argv[2] 

	#Remove log file, if exists
	if os.path.isfile(log_file_structures_not_accepted_by_pdb2gmx):
		os.remove(log_file_structures_not_accepted_by_pdb2gmx)

	#Get pdb files in ALL subdirectory
	my_pdb_files = get_files_pdb(path_pdb)

	for pdb_file in my_pdb_files:
		print pdb_file
		rename_atoms_structue_by_pdb2gmx(pdb_file, gmx_path)

main()	