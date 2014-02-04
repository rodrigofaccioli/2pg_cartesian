#
# --------------------------------------------------
# Author - Rodrigo Antonio Faccioli
# Date   - 27/01/2014
#
# --------------------------------------------------


CXX=g++
OPT=  -c -g
OPT_GSL = -static -Wall -lgsl -lgslcblas -lm -pedantic 

 all:protpred-Gromacs-Mono_PSP
 
protpred-Gromacs-Mono_PSP:protpred-Gromacs-Mono.o LoadConfig.o ea_mono.o topology.o pdbio.o protein.o futil.o pdbatom.o messages.o topologyio.o topologylib.o ea.o randomlib.o vector_math.o string_owner.o math_owner.o  osutil.o random_number_gsl.o load_parameters.o objective.o aminoacids.o aminoacids_io.o populationio.o
	$(CXX) -o protpred-Gromacs-Mono_PSP protpred-Gromacs-Mono.o LoadConfig.o ea_mono.o topology.o pdbio.o protein.o futil.o pdbatom.o messages.o topologyio.o topologylib.o ea.o randomlib.o vector_math.o string_owner.o math_owner.o  osutil.o random_number_gsl.o load_parameters.o objective.o aminoacids.o aminoacids_io.o populationio.o $(OPT_GSL) 

protpred-Gromacs-Mono.o:	protpred-Gromacs-Mono.c
	$(CXX) $(OPT) protpred-Gromacs-Mono.c -I../include

populationio.o:	populationio.c
	$(CXX) $(OPT) populationio.c -I../include

LoadConfig.o:	LoadConfig.cpp
	$(CXX) $(OPT) LoadConfig.cpp -I../include

topology.o:	topology.c
	$(CXX) $(OPT) topology.c -I../include

topologyio.o:	topologyio.c
	$(CXX) $(OPT) topologyio.c -I../include

topologylib.o:	topologylib.c
	$(CXX) $(OPT) topologylib.c -I../include

pdbio.o:	pdbio.c
	$(CXX) $(OPT) pdbio.c -I../include

energy_gromacs.o:	energy_gromacs.c
	$(CXX) $(OPT) energy_gromacs.c -I../include

ea_mono.o:	ea_mono.c
	$(CXX) $(OPT) ea_mono.c -I../include

ea.o:	ea.c
	$(CXX) $(OPT) ea.c -I../include

protein.o:	protein.c
	$(CXX) $(OPT) protein.c -I../include
	
futil.o:	futil.c
	$(CXX) $(OPT) futil.c -I../include

pdbatom.o:	pdbatom.c
	$(CXX) $(OPT) pdbatom.c -I../include

messages.o:	messages.c
	$(CXX) $(OPT) messages.c -I../include

randomlib.o:	randomlib.c
	$(CXX) $(OPT) randomlib.c -I../include

vector_math.o:	vector_math.c
	$(CXX) $(OPT) vector_math.c -I../include

string_owner.o:	string_owner.c
	$(CXX) $(OPT) string_owner.c -I../include

math_owner.o:	math_owner.c
	$(CXX) $(OPT) math_owner.c -I../include
	
random_number_gsl.o:	random_number_gsl.c
	$(CXX) $(OPT) random_number_gsl.c -I../include

osutil.o:	osutil.c
	$(CXX) $(OPT) osutil.c -I../include

load_parameters.o:	load_parameters.c
	$(CXX) $(OPT) load_parameters.c -I../include
	
ea_multi_obj_lib.o:	ea_multi_obj_lib.c
	$(CXX) $(OPT) ea_multi_obj_lib.c -I../include	

objective.o:	objective.c
	$(CXX) $(OPT) objective.c -I../include	

objective_io.o:	objective_io.c
	$(CXX) $(OPT) objective_io.c -I../include	

aminoacids.o:	aminoacids.c
	$(CXX) $(OPT) aminoacids.c -I../include	

aminoacids_io.o:	aminoacids_io.c
	$(CXX) $(OPT) aminoacids_io.c -I../include	

