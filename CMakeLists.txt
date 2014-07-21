cmake_minimum_required (VERSION 2.8)

project (2pg_cartesian CXX)

# Compiler specific settings:
if ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
		"${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	add_definitions (-D__const=const)
	add_definitions (-Dsnprintf=_snprintf)
	add_definitions (-Dgetpid=_getpid)
	add_definitions (-Dgetcwd=_getcwd)
	add_definitions (-Disnan=_isnan)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _CRT_SECURE_NO_WARNINGS")
endif ()

include_directories (${2pg_cartesian_SOURCE_DIR}/include)

add_subdirectory (src)
