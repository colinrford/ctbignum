#
# This file is part of
#
# CTBignum
#
# C++ Library for Compile-Time and Run-Time Multi-Precision and Modular Arithmetic
#
#
# This file is distributed under the Apache License, Version 2.0. See the LICENSE
# file for details.

cmake_minimum_required(VERSION 3.31.6)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.3)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.0)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 3.31.6)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
else()
  message(FATAL_ERROR "Version lower than 3.31.6 not supported (by me)")
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)

##
## PROJECT
## name and version
##
project(lam_ctbignum VERSION 0.7.1 LANGUAGES CXX)

##
## OPTIONS
##
option(LAM_CTBIGNUM_BuildTests "Build the unit tests when BUILD_TESTING is enabled." OFF)
option(LAM_CTBIGNUM_BuildBenchmarks "Build the benchmarks." OFF)

##
## CONFIGURATION
##
set(LAM_CTBIGNUM_TARGET_NAME               ${PROJECT_NAME})
set(LAM_CTBIGNUM_INCLUDE_BUILD_DIR "${PROJECT_SOURCE_DIR}/include/")

##
## TARGET
## create target and add module sources
##
add_library(${LAM_CTBIGNUM_TARGET_NAME} STATIC)

target_sources(${LAM_CTBIGNUM_TARGET_NAME}
    PUBLIC
    FILE_SET CXX_MODULES
    TYPE CXX_MODULES
    BASE_DIRS ${PROJECT_SOURCE_DIR}/include
    FILES
        # Primary module interface
        include/ctbignum/ctbignum.cppm
        # Module partitions (in dependency order)
        include/ctbignum/bigint.cppm
        include/ctbignum/slicing.cppm
        include/ctbignum/utility.cppm
        include/ctbignum/type_traits.cppm
        include/ctbignum/relational.cppm
        include/ctbignum/addition.cppm
        include/ctbignum/bitshift.cppm
        include/ctbignum/mult.cppm
        include/ctbignum/division.cppm
        include/ctbignum/gcd.cppm
        include/ctbignum/mod_inv.cppm
        include/ctbignum/barrett.cppm
        include/ctbignum/invariant_div.cppm
        include/ctbignum/montgomery.cppm
        include/ctbignum/mod_exp.cppm
        include/ctbignum/io.cppm
        include/ctbignum/literals.cppm
        include/ctbignum/field.cppm
)

add_library(lam::ctbignum ALIAS ${LAM_CTBIGNUM_TARGET_NAME})

##
## TESTS
## create and configure the unit test target
##
if(LAM_CTBIGNUM_BuildTests)
    include(CTest)
    if(BUILD_TESTING)
        enable_testing()
        add_subdirectory(test)
    endif()
endif()

if(LAM_CTBIGNUM_BuildBenchmarks)
    add_subdirectory(benchmarks)
endif()

##
## INSTALL
## install targets
##
include(GNUInstallDirs)

install(TARGETS ${LAM_CTBIGNUM_TARGET_NAME}
    EXPORT ${PROJECT_NAME}Targets
    FILE_SET CXX_MODULES
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
