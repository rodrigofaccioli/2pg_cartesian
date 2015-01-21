""" 
    Routines to join all pdb files into one pdb file
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import sys


""" This function obtains all pdb files 
    in mypath 
"""
def get_PROT_IND_files_pdb(mypath):
	only_pdb_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			#if file.endswith(".pdb"):
			if file.find("PROT_IND_") >=0:
				f_path = os.path.join(root,file)
				only_pdb_file.append(f_path)			
	return only_pdb_file

def main():

	main_pdb_path = sys.argv[1]

	model = int(1)
	str_model = "MODEL        #### \n"

	all_pdbs_path_file_name = os.path.join(main_pdb_path, "all_pdb.pdb")
	all_pdbs_file = open(all_pdbs_path_file_name, "w")

	pdb_files = get_PROT_IND_files_pdb(main_pdb_path)

	for pdb in pdb_files:		
		#preparing Model
		aux_model = str_model.replace("####", str(model))		
		all_pdbs_file.write(aux_model)
		#Reading pdb files
		aux_pdb_file  = open(pdb, "r")
		all_pdbs_file.write(aux_pdb_file.read())
		aux_pdb_file.close()
		#Finishing Model
		all_pdbs_file.write("TER \n")
		all_pdbs_file.write("ENDMDL \n")
		#Next model
		model = model + 1
	all_pdbs_file.close()

main()
