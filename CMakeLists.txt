#-------------------------------------------------------------------------------
# SuiteSparse/CMakeLists.txt:  root CMake build rules
#-------------------------------------------------------------------------------

# Copyright (c) 2023, Timothy A. Davis, All Rights Reserved.
# Just this particular file is under the Apache-2.0 license; each package has
# its own license.
# SPDX-License-Identifier: Apache-2.0

# The actually required minimum CMake version might be higher depending on the
# selected SuiteSparse projects.
cmake_minimum_required ( VERSION 3.20 )

project ( "SuiteSparse" )

#-------------------------------------------------------------------------------
# SuiteSparse policies
#-------------------------------------------------------------------------------

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/SuiteSparse_config/cmake_modules )

option ( ENABLE_CUDA "Enable CUDA acceleration" ON )

# SuiteSparsePolicy enables SUITESPARSE_CUDA if a CUDA compiler can be found.
include ( SuiteSparsePolicy )

#-------------------------------------------------------------------------------
# build options
#-------------------------------------------------------------------------------

# lower-case list of all projects that can be built by this root CMake file
set ( SUITESPARSE_ALL_PROJECTS
    "suitesparse_config;graphblas;lagraph" )

# lower-case list of extra projects that can be built by this root CMake file
set ( SUITESPARSE_EXTRA_PROJECTS
    "csparse" )

# lower-case list of known projects that can be built by this root CMake file
set ( SUITESPARSE_KNOWN_PROJECTS "${SUITESPARSE_ALL_PROJECTS};${SUITESPARSE_EXTRA_PROJECTS}" )

set ( SUITESPARSE_ENABLE_PROJECTS "all" CACHE STRING
    "Semicolon-separated list of SuiteSparse projects to be built (${SUITESPARSE_KNOWN_PROJECTS}, or \"all\")" )

# expand "all" early on
if ( "all" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    set ( SUITESPARSE_ENABLE_PROJECTS "${SUITESPARSE_ENABLE_PROJECTS};${SUITESPARSE_ALL_PROJECTS}" )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "all" )
    list ( REMOVE_DUPLICATES  SUITESPARSE_ENABLE_PROJECTS )
endif ( )

# check for unknown projects in list
foreach ( proj ${SUITESPARSE_ENABLE_PROJECTS} )
  if ( NOT "${proj}" IN_LIST SUITESPARSE_KNOWN_PROJECTS )
     message ( FATAL_ERROR "${proj} is not a known project: ${SUITESPARSE_KNOWN_PROJECTS}." )
  endif ( )
endforeach ( )

# Use OpenMP in CHOLMOD, GraphBLAS, LAGraph, and SuiteSparse_config
option ( NOPENMP "ON: do not use OpenMP.  OFF (default): use OpenMP" OFF )

# CHOLMOD options affecting dependencies
option ( NCAMD "ON: do not use CAMD/CCOLAMD.  OFF (default): use CAMD/CCOLAMD" OFF )

# KLU and UMFPACK options affecting dependencies
option ( NCHOLMOD "ON: do not use CHOLMOD in KLU and UMFPACK.  OFF (default): use CHOLMOD" OFF )

# overwrite BUILD_STATIC_LIBS specifically for GraphBLAS because building the
# library takes a long time
option ( GRAPHBLAS_BUILD_STATIC_LIBS "OFF (default): Do not build static libraries for GraphBLAS project.  ON: Use same value of BUILD_STATIC_LIBS for GraphBLAS like in the other projects" OFF )

# options to build with libraries installed on the system instead of building
# dependencies automatically
option ( USE_SYSTEM_GRAPHBLAS "ON: use GraphBLAS libraries installed on the build system.  OFF (default): Automatically build GraphBLAS as dependency if needed." OFF )
option ( USE_SYSTEM_SUITESPARSE_CONFIG "ON: use SuiteSparse_config libraries installed on the build system.  OFF (default): Automatically build SuiteSparse_config as dependency if needed." OFF )

#-------------------------------------------------------------------------------
# global variables
#-------------------------------------------------------------------------------

# Set to indicate that we are building from a root CMake file.
# That will change the directory layout and (imported) target names (namespace!)
# during the build process.
set ( SUITESPARSE_ROOT_CMAKELISTS ON )

#-------------------------------------------------------------------------------
# common SuiteSparse modules
#-------------------------------------------------------------------------------

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/SuiteSparse_config/cmake_modules )

include ( SuiteSparsePolicy )

if ( NOPENMP )
    # OpenMP has been disabled
    message ( STATUS "OpenMP disabled" )
    set ( OpenMP_C_FOUND OFF )
else ( )
    find_package ( OpenMP )
endif ( )

#-------------------------------------------------------------------------------
# check/add project dependencies
#-------------------------------------------------------------------------------

if ( USE_SYSTEM_GRAPHBLAS )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "graphblas" )
    find_package ( GraphBLAS REQUIRED )
else ( )
    if ( "lagraph" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT "graphblas" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"graphblas\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "graphblas" )
        endif ( )
    endif ( )
endif ( )

if ( USE_SYSTEM_BTF )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "btf" )
    find_package ( BTF REQUIRED )
else ( )
    if ( "klu" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT "btf" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"btf\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "btf" )
        endif ( )
    endif ( )
endif ( )

if ( USE_SYSTEM_CHOLMOD )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "cholmod" )
    find_package ( CHOLMOD REQUIRED )
else ( )
    if ( ( NOT NCHOLMOD AND
              ( "klu" IN_LIST SUITESPARSE_ENABLE_PROJECTS
                OR "umfpack" IN_LIST SUITESPARSE_ENABLE_PROJECTS ) )
            OR "spqr" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "paru" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT "cholmod" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"cholmod\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "cholmod" )
        endif ( )
    endif ( )
endif ( )

if ( USE_SYSTEM_AMD )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "amd" )
    find_package ( AMD REQUIRED )
else ( )
    if ( "cholmod" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "ldl" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "umfpack" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "spex" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT "amd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"amd\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "amd" )
        endif ( )
    endif ( )
endif ( )

if ( USE_SYSTEM_COLAMD )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "colamd" )
    find_package ( COLAMD REQUIRED )
else ( )
    if ( "cholmod" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "spex" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT "colamd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"colamd\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "colamd" )
        endif ( )
    endif ( )
endif ( )

if ( USE_SYSTEM_CAMD )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "camd" )
    find_package ( CAMD REQUIRED )
else ( )
    if ( NOT NCAMD AND "cholmod" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT USE_SYSTEM_CAMD AND NOT "camd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"camd\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "camd" )
        endif ( )
    endif ( )
endif ( )

if ( USE_SYSTEM_CCOLAMD )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "ccolamd" )
    find_package ( CCOLAMD REQUIRED )
else ( )
    if ( NOT NCAMD AND "cholmod" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT "ccolamd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"ccolamd\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "ccolamd" )
        endif ( )
    endif ( )
endif ( )

if ( USE_SYSTEM_SUITESPARSE_CONFIG )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "suitesparse_config" )
    find_package ( SuiteSparse_config REQUIRED )
else ( )
    if ( "mongoose" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "amd" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "btf" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "camd" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "ccolamd" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "colamd" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "cholmod" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "cxsparse" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "ldl" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "klu" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "umfpack" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "paru" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "rbio" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "spqr" IN_LIST SUITESPARSE_ENABLE_PROJECTS
            OR "spex" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        if ( NOT "suitesparse_config" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
            message ( STATUS "Adding \"suitesparse_config\" to the list of built targets." )
            list ( APPEND SUITESPARSE_ENABLE_PROJECTS "suitesparse_config" )
        endif ( )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# include selected projects
#-------------------------------------------------------------------------------

if ( "suitesparse_config" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    add_subdirectory ( "SuiteSparse_config" )
    if ( TARGET SuiteSparseConfig )
        add_library ( SuiteSparse::SuiteSparseConfig ALIAS SuiteSparseConfig )
    else ( )
        add_library ( SuiteSparse::SuiteSparseConfig ALIAS SuiteSparseConfig_static )
    endif ( )
    if ( TARGET SuiteSparseConfig_static )
        add_library ( SuiteSparse::SuiteSparseConfig_static ALIAS SuiteSparseConfig_static )
    endif ( )
endif ( )

if ( "graphblas" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    add_subdirectory ( "GraphBLAS" )
    if ( TARGET GraphBLAS )
        add_library ( SuiteSparse::GraphBLAS ALIAS GraphBLAS )
    else ( )
        add_library ( SuiteSparse::GraphBLAS ALIAS GraphBLAS_static )
    endif ( )
    if ( TARGET GraphBLAS_static )
        add_library ( SuiteSparse::GraphBLAS_static ALIAS GraphBLAS_static )
    endif ( )
endif ( )

if ( "lagraph" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    add_subdirectory ( "LAGraph" )
    if ( TARGET LAGraph )
        add_library ( SuiteSparse::LAGraph ALIAS LAGraph )
    else ( )
        add_library ( SuiteSparse::LAGraph ALIAS LAGraph_static )
    endif ( )
    if ( TARGET LAGraph_static )
        add_library ( SuiteSparse::LAGraph_static ALIAS LAGraph_static )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# report status
#-------------------------------------------------------------------------------

include ( SuiteSparseReport )

#-------------------------------------------------------------------------------
# enable testing facilities
#-------------------------------------------------------------------------------

# Currently, only LAGraph, Mongoose, and CHOLMOD have ctests.

# FIXME: convert more of the existing demos to ctests.

# Most packages have a ./Tcov folder with a full statement coverage test,
# but these are not imported into cmake yet.

# Most packages also have a ./Demo folder, with shorter examples.  These would
# be nice to add as quick ctests.

# CHOLMOD/Tcov takes about 20 minutes to run.  It is also a full coverage
# test of AMD, CAMD, COLAMD, and CCOLAMD, however.  The current CHOLMOD
# ctest is based on a few ./Demo programs.  It's fast but not a full coverate
# test.

# The CSparse/CXSparse Tcov tests are very fast and would be good candidates to
# add.

enable_testing ( )

