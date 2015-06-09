cmake_minimum_required(VERSION 2.8)

# project Information 
project(2pg_cartesian)
set(PROJECT_VERSION "1.0")


# Compiler specific settings:
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	add_definitions (-D__const=const)
	add_definitions (-Dsnprintf=_snprintf)
	add_definitions (-Dgetpid=_getpid)
	add_definitions (-Dgetcwd=_getcwd)
	add_definitions (-Disnan=_isnan)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _CRT_SECURE_NO_WARNINGS")
endif ()

# Set compiler flags
SET ( CMAKE_CXX_FLAGS "-lm -pedantic")

#Set CXX compiler for all files below
set_source_files_properties(include/LoadConfig.h PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-NSGA2.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Dominance.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Sort_Method_by_Front_Dominance.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-MC_Metropolis.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Mono.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Random_Algorithm.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Sort_Method_Files_by_Front_Dominance.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/LoadConfig.cpp  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/ea_mono.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/topology.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/pdbio.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protein.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/futil.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/pdbatom.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/messages.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/topologyio.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/topologylib.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/randomlib.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/vector_math.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/string_owner.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/math_owner.c   PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/osutil.c   PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/load_parameters.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/objective.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/aminoacids.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/aminoacids_io.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/populationio.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/rotation.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/solution.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/gromacs_objectives.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/solutionio.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/algorithms.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/ea_nsga2.c  PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/dominance.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/owner_file_analysis.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/mc_metropolis.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/random_algorithm.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/diehdral.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Front.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Crossover.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Test_compute_Diehdral.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/protpred-Gromacs-Test_rotation.c PROPERTIES LANGUAGE CXX )

# set include
include_directories(include)

# add libries
add_library(2PG_lib STATIC 
src/LoadConfig.cpp 
src/ea_mono.c 
src/topology.c 
src/pdbio.c 
src/protein.c 
src/futil.c 
src/pdbatom.c 
src/messages.c 
src/topologyio.c 
src/topologylib.c 
src/randomlib.c 
src/vector_math.c 
src/string_owner.c 
src/math_owner.c  
src/osutil.c  
src/load_parameters.c 
src/objective.c 
src/aminoacids.c 
src/aminoacids_io.c 
src/populationio.c 	
src/rotation.c 
src/solution.c 
src/gromacs_objectives.c 
src/solutionio.c 
src/algorithms.c 
src/ea_nsga2.c 
src/dominance.c
src/owner_file_analysis.c
src/mc_metropolis.c
src/random_algorithm.c
src/diehdral.c
) #end of 2PG_lib

# add target
add_executable(protpred-Gromacs-NSGA2 src/protpred-Gromacs-NSGA2.c)
target_link_libraries(protpred-Gromacs-NSGA2 2PG_lib)

add_executable(protpred-Gromacs-Dominance src/protpred-Gromacs-Dominance.c)
target_link_libraries(protpred-Gromacs-Dominance 2PG_lib)

add_executable(protpred-Gromacs-Front src/protpred-Gromacs-Front.c)
target_link_libraries(protpred-Gromacs-Front 2PG_lib)

add_executable(protpred-Gromacs-Sort_Method_by_Front_Dominance src/protpred-Gromacs-Sort_Method_by_Front_Dominance.c)
target_link_libraries(protpred-Gromacs-Sort_Method_by_Front_Dominance 2PG_lib)

add_executable(protpred-Gromacs-MC_Metropolis src/protpred-Gromacs-MC_Metropolis.c)
target_link_libraries(protpred-Gromacs-MC_Metropolis 2PG_lib)

add_executable(protpred-Gromacs-Mono src/protpred-Gromacs-Mono.c)
target_link_libraries(protpred-Gromacs-Mono 2PG_lib)

add_executable(protpred-Gromacs-Random_Algorithm src/protpred-Gromacs-Random_Algorithm.c)
target_link_libraries(protpred-Gromacs-Random_Algorithm 2PG_lib)

add_executable(protpred-Gromacs-Sort_Method_Files_by_Front_Dominance src/protpred-Gromacs-Sort_Method_Files_by_Front_Dominance.c)
target_link_libraries(protpred-Gromacs-Sort_Method_Files_by_Front_Dominance 2PG_lib)

add_executable(protpred-Gromacs-Crossover src/protpred-Gromacs-Crossover.c)
target_link_libraries(protpred-Gromacs-Crossover 2PG_lib)

add_executable(protpred-Gromacs-Test_compute_Diehdral src/protpred-Gromacs-Test_compute_Diehdral.c)
target_link_libraries(protpred-Gromacs-Test_compute_Diehdral 2PG_lib)

add_executable(protpred-Gromacs-Test_rotation src/protpred-Gromacs-Test_rotation.c)
target_link_libraries(protpred-Gromacs-Test_rotation 2PG_lib)


# install
install(TARGETS protpred-Gromacs-NSGA2 DESTINATION bin)
install(TARGETS protpred-Gromacs-Dominance DESTINATION bin)
install(TARGETS protpred-Gromacs-Front DESTINATION bin)
install(TARGETS protpred-Gromacs-MC_Metropolis DESTINATION bin)
install(TARGETS protpred-Gromacs-Mono DESTINATION bin)
install(TARGETS protpred-Gromacs-Random_Algorithm DESTINATION bin)
install(TARGETS protpred-Gromacs-Sort_Method_Files_by_Front_Dominance DESTINATION bin)
install(TARGETS protpred-Gromacs-Sort_Method_by_Front_Dominance DESTINATION bin)
install(TARGETS protpred-Gromacs-Crossover DESTINATION bin)
install(TARGETS protpred-Gromacs-Test_compute_Diehdral DESTINATION bin)
install(TARGETS protpred-Gromacs-Test_rotation DESTINATION bin)

#Coping script directory
install(DIRECTORY scripts/ DESTINATION scripts 
	PATTERN "scripts/*" 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ
)#End of Coping script directory

#Coping test directory
install(DIRECTORY test/ DESTINATION test 
	PATTERN "test/*" 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ
)#End of Coping test directory