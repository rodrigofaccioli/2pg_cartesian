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
	algorithms.c)

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
