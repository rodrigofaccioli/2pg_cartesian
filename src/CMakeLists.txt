file (GLOB BASE_SOURCES
	LoadConfig.cpp
	topology.c
	pdbio.c
	protein.c
	futil.c
	pdbatom.c
	messages.c
	topologyio.c
	topologylib.c
	randomlib.c
	vector_math.c
	string_owner.c
	math_owner.c
	osutil.c
	load_parameters.c
	objective.c
	aminoacids.c
	aminoacids_io.c
	populationio.c
	rotation.c
	solution.c
#	gromacs_objectives.c
	solutionio.c
	algorithms.c)

set_source_files_properties (${BASE_SOURCES} PROPERTIES LANGUAGE CXX)

#set_source_files_properties (
#	protpred-Gromacs-Random_Algorithm.c
#	random_algorithm.c
#	PROPERTIES LANGUAGE CXX)

add_library (protpred-Gromacs-Random_Algorithm_PSP
	SHARED
	protpred-Gromacs-Random_Algorithm.c
	random_algorithm.c
	${BASE_SOURCES})

