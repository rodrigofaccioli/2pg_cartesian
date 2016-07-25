#! /usr/bin/env python

#This script extracts all front files
# python extract_front_files.py <directory where are front files>
# Example: python count_global_non_dominated_solution.py /home/faccioli/Execute/T0773/


import os
import sys
from shutil import copy

def get_files_front(mypath):
#Obtaing the tar.gz files from mypath
	only_tar_file = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file.endswith(".front"):
				f_path = os.path.join(root,file)
				only_tar_file.append(f_path)			
	return only_tar_file

def get_all_directories(mypath):
	only_directory = []
	for root, dirs, files in os.walk(mypath):
		for directory in dirs:
			only_directory.append(directory)
	return only_directory

def save_non_dominated(frontfile, result_front):
	obj1 = ""
	obj2 = ""
	obj_comb = ""
	f_front = open(frontfile)
	f_res_front = open(result_front,"a")
	for line in f_front.readlines():
		if line.find("#") == -1: # NO found character # means contain data
			front     = int(line.split()[1])
			# when is no-dominated
			if front == 0:
				dominated = line.split()[2]
				method    = line.split()[5]				
				outline = method + "\t\t" + str(dominated) + "\t" + obj_comb +"\n"
				f_res_front.write(outline)	
		elif line.find("#Ranking  Front") == 0:			
			obj1 = line.split()[3]
			obj2 = line.split()[4]
			obj_comb = obj1 + " x " + obj2

			
			
	f_res_front.close()
	f_front.close()


#Represents how many times method was non-dominated considering all combination of objectives
def count_global_non_dominated_method(result_front, count_front_method):
	dic = {}
	list_method = []
	f_front = open(result_front)
	for line in f_front.readlines():
		method = line.split()[0]
		list_method.append(method)
	f_front.close()

	f_count_front_method = open(count_front_method,"w")
	outline = "#Method" +"\t" + "Count"+"\n"
	f_count_front_method.write(outline)
	for count, elem in sorted(((list_method.count(e), e) for e in set(list_method)), reverse=True):
		outline = elem +"\t" + str(count)+"\n"
		f_count_front_method.write(outline) 
	f_count_front_method.close()

def main():

	path_root = sys.argv[1] #main path

	result_front = os.path.join(path_root, "result_front.txt")	
	if os.path.isfile(result_front):
		os.remove(result_front)

	count_front_method  = os.path.join(path_root, "result_count_global_non_dominated_method.txt")	

	#Get front file in all subdirectory
	my_front_files = get_files_front(path_root)

	#Copy front files to
	for front_file in my_front_files:
		save_non_dominated(front_file, result_front)

	count_global_non_dominated_method(result_front, count_front_method)

main()