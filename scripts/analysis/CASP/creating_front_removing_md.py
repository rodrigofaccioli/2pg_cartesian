#! /usr/bin/env python
""" 
    Routines to create front file removing md.1.last solutions.
    So, in this front file contains solution from CASP11 stage2
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""
#Example: creating_front_removing_md.py /home/faccioli/Dropbox/Rodrigo_FCFRP/CASP11/analises/2PG /home/faccioli/workspace/2pg_cartesian/build/protpred-Gromacs-Sort_Method_by_Front_Dominance

import sys
import os

""" This function obtains all sub-directories 
	from root_path
"""
def get_all_directory(root_path):
	list_res = []
	aux_direct = os.walk(root_path)
	for l in aux_direct:
		list_res.append(l[0])
	return list_res

""" This function obtains all front files 
    in mypath 
"""
def get_files_front(mypath):
	only_pdbqt_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".front"):
				f_path = os.path.join(root,file)
				only_pdbqt_file.append(f_path)			
	return only_pdbqt_file

def main():
	#root path
	root_path = sys.argv[1]
	#path_file 2PG program
	program_2PG = sys.argv[2]


	str_to_remove = ".md.1_last"
	input_file_sort = "front_without_md.txt"
	header_input_file_sort = ";-1	-1	0 \n"
	command_ref = "@PROGRAM@ @INPUT_FILE@ 2>/dev/null"

	all_directories = get_all_directory(root_path)
	for directory in all_directories:
		front_file = get_files_front(directory)
		#checking list is empty
		if len(front_file) > 0: #not empty
			#preparing filename to be saved
			full_input_file_sort = os.path.join(directory, input_file_sort)	
			input_file = open(full_input_file_sort, "w")
			#preparing filename to be read
			front_file_str = str(front_file[0])
			front_file = open(front_file_str, "r")
			#reading front file
			input_file.write(header_input_file_sort)
			for line in front_file.readlines():
				if str(line).find("#") < 0 and str(line).find(str_to_remove) < 0:
					line_splited = line.split()
					line_save = str(float(line_splited[3]))  +"\t"+ str(float(line_splited[4])) +"\t"+ str(line_splited[5]) +"\n"
					input_file.write(line_save)
			input_file.close()
			front_file.close()
			#Running 2PG program			
			os.chdir(directory)
			command = str(command_ref).replace("@PROGRAM@", program_2PG).replace("@INPUT_FILE@", full_input_file_sort)
			os.system(command)
main()