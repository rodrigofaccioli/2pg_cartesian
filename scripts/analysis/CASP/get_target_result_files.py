#! /usr/bin/env python
""" 
    Routines to download the CASP11 results from predictioncenter.org
    These routines were developed by:
    Rodrigo Antonio Faccioli - rodrigo.faccioli@usp.br / rodrigo.faccioli@gmail.com  
    Leandro Oliveira Bortot  - leandro.bortot@usp.br / leandro.obt@gmail.com 
"""

import sys
import urllib
import os

#example get_target_result_files.py <main_path> <all_targets_file>
#    get_target_result_files.py /home/faccioli/Dropbox/Rodrigo_FCFRP/CASP11/analises/2PG/ /home/faccioli/Dropbox/Rodrigo_FCFRP/CASP11/analises/2PG/2pg_casp11_all_targets.txt

def get_target_name(a):
	s1 = str(a).split("TS")[0]

	if str(a).find("D") >=0:
		s2 = str(a).split("-")[1]
		target = s1+"-"+s2
	else:	
		target = s1
	return str(target).strip()

def main():

	ref_dir = str(sys.argv[1])
	all_targets_file = str(sys.argv[2])
	#main url
	url_main="http://www.predictioncenter.org/casp11/results.cgi?dm_class=all&model=&pc=&access_type=1&multi_sort=&groups_id=&target=#####&targets_list=&result_id=&tr_type=all&order=&groups_list=&results=all&view=txt"		

	file_tar = open(all_targets_file, "r")
	for line in file_tar:
		target = get_target_name(line)
		#creating directory
		full_path = os.path.join(ref_dir,target)
		os.mkdir(full_path)

		#replace string ##### for target
		url=str(url_main).replace("#####",target)

		#retrieving file and saving it as casp11.txt
		retrieving_file = os.path.join(full_path,"casp11.txt")
		urllib.urlretrieve (url, retrieving_file)


main()