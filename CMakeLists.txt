cmake_minimum_required (VERSION 2.8)

project (2pg_cartesian CXX)

# Compiler specific settings:
if ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
		"${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

endif ()

include_directories (${2pg_cartesian_SOURCE_DIR}/include)

add_subdirectory (src)
