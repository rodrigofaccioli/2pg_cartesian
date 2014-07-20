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
#	gromacs_objectives.c
	solutionio.c
	algorithms.c)

set_source_files_properties (${BASE_SOURCES} PROPERTIES LANGUAGE CXX)

add_library (2pg_cartesian
	SHARED
	${BASE_SOURCES})


# Declare macros to export DLL functions:
include (GenerateExportHeader)

GENERATE_EXPORT_HEADER (
	2pg_cartesian
	BASE_NAME 2pg_cartesian)

include_directories (${CMAKE_CURRENT_BINARY_DIR})

# This is in order to export rather than import things tagged with
# _2PG_CARTESIAN_EXPORT
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	add_definitions (-D_2pg_cartesian_EXPORTS)
endif ()


set_source_files_properties (
	protpred-Gromacs-Random_Algorithm.c random_algorithm.c
	PROPERTIES LANGUAGE CXX)

add_executable (protpred-Gromacs-Random_Algorithm_PSP
	protpred-Gromacs-Random_Algorithm.c random_algorithm.c)

target_link_libraries (protpred-Gromacs-Random_Algorithm_PSP
	2pg_cartesian)
