""" 
    Routines to compute RMSD of all PROT_IND_ files
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import os
import sys
from collections import OrderedDict

native = "1VII.pdb"
path_gromacs ="/home/faccioli/Programs/gmx-4.6.5/no_mpi/bin/"


main_command = "echo C-alpha C-alpha | @PATH_GROMACS@./g_rms -f @PROT@ -s @NATIVE@ -o temporary_rmsd.xvg 2>/dev/null"

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
	pdb_path = sys.argv[1]

	dict_rmsd = {}
	all_pdbs = get_PROT_IND_files_pdb(pdb_path)
	for pdb in all_pdbs:
		aux_command = main_command.replace("@PATH_GROMACS@", path_gromacs).replace("@PROT@",pdb).replace("@NATIVE@", native)
		os.system(aux_command)
		temp_rmsd = open("temporary_rmsd.xvg", "r")
		for line in temp_rmsd.readlines():
			if line.find("@") < 0 and line.find("#") <0:				
				rmsd_value = float(str(line).split()[1])
				only_pdb_file_name = os.path.basename(pdb)
				dict_rmsd[only_pdb_file_name] = rmsd_value
		temp_rmsd.close()
		os.remove("temporary_rmsd.xvg")
	#Saving dictionary
	rmsd_final = open("all_rmsd.txt", "w")
	d_sorted_by_value = OrderedDict(sorted(dict_rmsd.items(), key=lambda x: x[1]))
	for key, value in d_sorted_by_value.items():
		rmsd_final.write(str(key) +"\t" + str(value) + "\n")
	rmsd_final.close()

main()