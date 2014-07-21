###############################
## - Create library target - ##
###############################

# shared library sources
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
	gromacs_objectives.c
	solutionio.c
	algorithms.c
	ea_mono.c
	ea_nsga2.c
	dominance.c
#	owner_file_analysis.c # Missing dirent.h opendir etc implementation for Windows
	)

# This is in order to export rather than import things tagged with
# _2PG_CARTESIAN_EXPORT.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	add_definitions (-D_2pg_cartesian_EXPORTS)
endif ()

# Specify that the sources are inteded to be compiled as C++
set_source_files_properties (${BASE_SOURCES} PROPERTIES LANGUAGE CXX)

# Add library 2pg_cartesian target
add_library (2pg_cartesian
	SHARED
	${BASE_SOURCES})


# Declare macros to export DLL functions:
include (GenerateExportHeader)

# All exported functions anc classes must be tagged with the
# 2PG_CARTESIAN_EXPORT macro defined in the genearated 2pg_cartesian_export.h
# header file:
GENERATE_EXPORT_HEADER (
	2pg_cartesian
	BASE_NAME 2pg_cartesian)

# The current binary directory must be in the path to look for header files, as
# 2pg_cartesian_export.h is generated here:
include_directories (${CMAKE_CURRENT_BINARY_DIR})



###################################
## - Create executable targets - ##
###################################

set_source_files_properties (
	protpred-Gromacs-Random_Algorithm.c random_algorithm.c
	PROPERTIES LANGUAGE CXX)

add_executable (protpred-Gromacs-Random_Algorithm_PSP
	protpred-Gromacs-Random_Algorithm.c random_algorithm.c)

target_link_libraries (protpred-Gromacs-Random_Algorithm_PSP
	2pg_cartesian)



set_source_files_properties (
	protpred-Gromacs-MC_Metropolis.c mc_metropolis.c
	PROPERTIES LANGUAGE CXX)

add_executable (protpred-Gromacs-MC_Metropolis_PSP
	protpred-Gromacs-MC_Metropolis.c mc_metropolis.c)

target_link_libraries (protpred-Gromacs-MC_Metropolis_PSP
	2pg_cartesian)



set_source_files_properties (
	protpred-Gromacs-Mono.c
	PROPERTIES LANGUAGE CXX)

add_executable (protpred-Gromacs-Mono_PSP
	protpred-Gromacs-Mono.c)

target_link_libraries (protpred-Gromacs-Mono_PSP
	2pg_cartesian)



set_source_files_properties (
	protpred-Gromacs-Dominance.c
	PROPERTIES LANGUAGE CXX)

add_executable (protpred-Gromacs-Dominance
	protpred-Gromacs-Dominance.c)

target_link_libraries (protpred-Gromacs-Dominance
	2pg_cartesian)



set_source_files_properties (
	protpred-Gromacs-NSGA2.c
	PROPERTIES LANGUAGE CXX)

add_executable (protpred-Gromacs-NSGA2_PSP
	protpred-Gromacs-NSGA2.c)

target_link_libraries (protpred-Gromacs-NSGA2_PSP
	2pg_cartesian)



# The following targets are commented out because of missing functions from
# owner_file_analysis.c which is excluded from the 2pg_cartesian_library for
# lack of dirent.h, opendir etc. implementation on Windows:

#set_source_files_properties (
#	protpred-Gromacs-Front.c
#	PROPERTIES LANGUAGE CXX)

#add_executable (protpred-Gromacs-Front
#	protpred-Gromacs-Front.c)

#target_link_libraries (protpred-Gromacs-Front
#	2pg_cartesian)



#set_source_files_properties (
#	protpred-Gromacs-Sort_Method_Files_by_Front_Dominance.c
#	PROPERTIES LANGUAGE CXX)

#add_executable (protpred-Gromacs-Sort_Method_Files_by_Front_Dominance
#	protpred-Gromacs-Sort_Method_Files_by_Front_Dominance.c)

#target_link_libraries (protpred-Gromacs-Sort_Method_Files_by_Front_Dominance
#	2pg_cartesian)
