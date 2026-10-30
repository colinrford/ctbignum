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

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LamCxxModuleStd.cmake)

## PROJECT
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" _lam_ctbignum_version)
string(STRIP "${_lam_ctbignum_version}" _lam_ctbignum_version)
project(lam_ctbignum VERSION ${_lam_ctbignum_version} LANGUAGES CXX)

set_property(DIRECTORY APPEND PROPERTY
  CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LamModule.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Sanitizers.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Warnings.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LamConstexprLimit.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CtbignumHelpers.cmake)

## OPTIONS
# The Conan recipe turns both off; test/ and benchmarks/ aren't exported.
option(LAM_CTBIGNUM_BUILD_TESTS "Build the unit tests." ON)
option(LAM_CTBIGNUM_BUILD_BENCHMARKS "Build the benchmarks." ON)

## TARGET
add_library(lam_ctbignum STATIC)
target_sources(lam_ctbignum
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
        include/ctbignum/decimal_literals.cppm
        include/ctbignum/field.cppm
        include/ctbignum/roots.cppm
)
add_library(lam::ctbignum ALIAS lam_ctbignum)
add_library(lam_ctbignum::ctbignum ALIAS lam_ctbignum)

# Target stays lam_ctbignum: ctbignum_add_executable hard-codes lam_ctbignum.dir.
set_target_properties(lam_ctbignum PROPERTIES
  OUTPUT_NAME lam_ctbignum
  EXPORT_NAME ctbignum)

target_compile_features(lam_ctbignum PUBLIC cxx_std_23)
lam_target_warnings(lam_ctbignum)

# Not check_cxx_source_compiles: try-compile inherits CMAKE_CXX_MODULE_STD and fails.
if(MSVC)
  set(CBN_HAS_INT128_BOOL "false")
  set(CBN_UINT128_TYPE_ALIAS "using uint128_or_void = void;")
else()
  set(CBN_HAS_INT128_BOOL "true")
  set(CBN_UINT128_TYPE_ALIAS "__extension__ using uint128_or_void = unsigned __int128;")
endif()

lam_add_config_partition(lam_ctbignum NAME ctbignum
  TEMPLATE ${CMAKE_CURRENT_SOURCE_DIR}/include/ctbignum/ctbignum_config.cppm.in
  EXTRA_DEFINES
    -DCBN_HAS_INT128_BOOL=${CBN_HAS_INT128_BOOL}
    "-DCBN_UINT128_TYPE_ALIAS=${CBN_UINT128_TYPE_ALIAS}")

enable_testing()
lam_add_config_dump(lam_ctbignum NAME ctbignum)

## TESTS / BENCHMARKS
if(LAM_CTBIGNUM_BUILD_TESTS)
    add_subdirectory(test)
endif()

if(LAM_CTBIGNUM_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

## INSTALL
lam_install_module_package(lam_ctbignum
  NAME        ctbignum
  DESCRIPTION "lam.ctbignum — compile-time big integer arithmetic (C++23 modules)"
  LICENSE     "Apache-2.0"
)
