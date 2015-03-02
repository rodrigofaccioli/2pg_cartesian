#! /usr/bin/env python
""" 
    Routines to obtain all ranks of branch for each target
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

#Example:  get_all_ranks_of_group_name_each_target.py /home/faccioli/Execute/2PG/teste_casp11/ 2PG casp11.txt

import os
import sys
from collections import OrderedDict

""" This function obtains all sub-directories 
	from root_path
"""
def get_all_directory(root_path):
	list_res = []
	aux_direct = os.walk(root_path)
	for l in aux_direct:
		list_res.append(l[0])
	return list_res

""" This function obtains file of casp results
"""
def get_file_casp_results(mypath, file_casp_results):
	only_file_casp_results = []
	for root, dirs, files in os.walk(mypath):
		for file in files:
			if file == file_casp_results:
				f_path = os.path.join(root,file)				
				only_file_casp_results.append(f_path)			
	return only_file_casp_results

""" This function returns a list of models and its ranks of group name
based on target
"""
def get_models_and_ranks_of_branch(path_casp_results, group_name):
	return_dic = {}	
	splited_line = []
	f_casp_results = open(path_casp_results)
	for line in f_casp_results.readlines():
		splited_line = str(line).split()
		if splited_line:
			if str(splited_line[3]) == group_name:
				model = str(splited_line[1])
				rank = int(splited_line[6])
				return_dic[model] = rank

	f_casp_results.close()

	return return_dic

""" This function creates the group name rank models 
"""
def create_file_group_name_rank_models(d_sorted_by_value, group_name, root_path):
	file_name_output = group_name+"_rank_models.txt" 
	header = ";Model        \t Rank"
	path_file_name_output = os.path.join(root_path, file_name_output)
	output_final = open(path_file_name_output, "w")
	output_final.write(header + "\n")
	for key, value in d_sorted_by_value.items():
		if value > 0:		
			output_final.write(str(key) +"      \t " + str(value) + "\n")
	output_final.close()


""" This function creates the histogram of group name rank models 
"""
def create_file_histogram_group_name_rank_models(d_sorted_by_value, group_name, root_path):

	#Building rank and frequency
	dict_file = {}
	ref_rank = int(d_sorted_by_value.values()[0])	
	dict_file[ref_rank] = 0
	for l_item in d_sorted_by_value.items():
		if int(l_item[1]) == ref_rank:			
				dict_file[ref_rank] = dict_file[ref_rank] + 1
		else:
			ref_rank = int(l_item[1])
			dict_file[ref_rank] = 1

	#Creating histogram file
	file_name_output = group_name+"_histogram_rank_models.xvg"
	header = "#Rank\tFrequency"
	path_file_name_output = os.path.join(root_path, file_name_output)
	output_final = open(path_file_name_output, "w")
	output_final.write(header+"\n")
	for key, value in dict_file.items():
		if key > 0:
			line = str(key)+"\t"+str(value)+"\n"		
			output_final.write(line)
	output_final.close()

def main():

	#root directory
	root_path = str(sys.argv[1])
	#branch name that want to know its ranks of each target
	group_name = str(sys.argv[2])
	#file where the results are. This file is created by get_target_result_files.py
	file_casp_results = str(sys.argv[3])

	dict_result = {}

	all_directories = get_all_directory(root_path)
	for directory in all_directories:
		path_casp_results = str(get_file_casp_results(directory, file_casp_results)[0])		
		dict_temp = get_models_and_ranks_of_branch(path_casp_results, group_name)
		if dict_temp:
			#means there is model of group name
			for key, value in  dict_temp.items():
				dict_result[key] = value

	#sorting dictionary for output files
	d_sorted_by_value = OrderedDict(sorted(dict_result.items(), key=lambda x: x[1]))

	#creating the output files
	create_file_group_name_rank_models(d_sorted_by_value, group_name, root_path)
	create_file_histogram_group_name_rank_models(d_sorted_by_value, group_name, root_path)

main()