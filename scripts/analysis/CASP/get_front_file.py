#! /usr/bin/env python
""" 
    Routines to obtain the front files that 
    were output of protpred-Gromacs-Sort_Method_Files_by_Front_Dominance program
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

#Example:  get_front_file.py /home/faccioli/Dropbox/Rodrigo_FCFRP/CASP11/analises/2PG/2pg_casp11_all_targets.txt /home/faccioli/Dropbox/Rodrigo_FCFRP/CASP11/analises/2PG/

import os
import sys
import tarfile

def extract_front_file_from_tar_gz(tar_gz_path_file, extract_dest_path):
	if os.path.isfile(tar_gz_path_file):
		tar = tarfile.open(tar_gz_path_file, 'r')
		for item in tar:
			if item.name.find('.front') >=0:			
				tar_file_name = os.path.basename(item.name)	
				path_file_name = os.path.join(extract_dest_path, tar_file_name)
				f = tar.extractfile(item.name)
				save_file = open(path_file_name,"w")
				save_file.write(f.read())
				save_file.close()
		tar.close()


def get_target_name(a):
	s1 = str(a).split("TS")[0]

	if str(a).find("D") >=0:
		s2 = str(a).split("-")[1]
		target = s1+"-"+s2
	else:	
		target = s1
	return str(target).strip()

def get_target_name_for_front_file(aux_target):
	if str(aux_target).find("D") >=0:
		target = str(aux_target).split("-")[0]		 
	else:	
		target = aux_target
	return str(target).strip()

def main():
	
	#File that contains all targets
	targets_file = sys.argv[1]
	#root path for saving front files
	extract_dest_path = sys.argv[2]
	#root path for searching front files
	search_front_main_path = sys.argv[3]

	ref_front_tar_gz="####_f_Pot_GBSA_1.tar.gz"
	
	all_targets_file = open(targets_file, "r")
	for target in all_targets_file:
		target_name = get_target_name(target)
		dest_full_path = os.path.join(extract_dest_path, target_name)
		#preparing the tar filename				
		target_name_for_obtain_front_filename = get_target_name_for_front_file(target_name)
		full_target_name_for_obtain_front_filename = os.path.join(search_front_main_path, target_name_for_obtain_front_filename)
		target_name_for_obtain_front_filename = target_name_for_obtain_front_filename[1:]		
		target_name_tar_gz = str(ref_front_tar_gz).replace("####",target_name_for_obtain_front_filename)
		tarfile = os.path.join(full_target_name_for_obtain_front_filename, target_name_tar_gz)
		#extracting front file		
		extract_front_file_from_tar_gz(tarfile, dest_full_path)		

main()