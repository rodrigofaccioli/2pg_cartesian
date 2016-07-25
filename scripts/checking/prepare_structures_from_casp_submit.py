#! /usr/bin/env python

# This script prepare structures of directory from CASP
# It checks name of atoms and termination of file
# IMPORTANT: Execute this script when model is read to be submitted to CASP.
# Therefore, full all fields of model (PFRMAT, TARGET, AUTHOR, REMARK, METHOD, MODEL and PARENT)
# BEFORE run this script.
#
# It uses bio.pdb module of BioPython project. 
# apt-get build-dep python-biopython is command to install lastest version of BioPython project
#
# prepare_structures_for_casp_submission.py <pdbs path> 
# Example: prepare_structures_from_casp_submission.py /home/faccioli/Dropbox/CASP11/T0773/


import sys
import os
from Bio.PDB import PDBParser, PDBIO

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

def checking_pdb_to_submit_casp(pdb_file):
	parser = PDBParser(PERMISSIVE=1)	
	structure = parser.get_structure('', pdb_file)

	for model in structure:
		for chain in model:	
			#Rename CD atom of ILE
			for residue in chain:
				if residue.resname == 'ILE':
					for atom in residue:
						if atom.get_name() == 'CD1':
							residue['CD1'].fullname = ' CD'						
	return structure

def saving_pdb(pdb_file, structure):
	w = PDBIO()
	w.set_structure(structure)
	w.save(pdb_file)


def get_header(pdb_file):
	header = []
	f_pdb = open(pdb_file)
	for line in f_pdb.readlines():
		if line.find("ATOM") < 0:
			header.append(line)
		else:
			break 
	f_pdb.close()
	return header

def creating_final_pdb(pdb_file, header):
	pdb_aux_name = "aux_pdb.pdb"
	f_pdb_aux = open(pdb_aux_name, "w")
	#Seting header
	for h in header:
		f_pdb_aux.write(h)
	#Coping file pdb
	f_pdb = open(pdb_file)
	for line in f_pdb.readlines():
		f_pdb_aux.write(line)
	f_pdb.close() #close file pdb
	#Adding END word in final of temp file
	f_pdb_aux.write("END")
	f_pdb_aux.close()

	#remove pdb final
	os.remove(pdb_file)

	#rename aux_pdb to pdb_file
	os.rename(pdb_aux_name, pdb_file)

def main():
	
	#path of pdb files
	path_pdb = sys.argv[1]

	#Get pdb files in ALL subdirectory
	my_pdb_files = get_files_pdb(path_pdb)

	for pdb_file in my_pdb_files:		
		header = get_header(pdb_file)				
		structure = checking_pdb_to_submit_casp(pdb_file)
		saving_pdb(pdb_file, structure)
		creating_final_pdb(pdb_file, header)

main()	
