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

#############################################################################
# Catch library with the main function to speed up build
#############################################################################

add_library(catch_main OBJECT
    "src/tests_main.cpp"
)
set_target_properties(catch_main PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED ON
)
target_include_directories(catch_main PRIVATE "thirdparty/catch")

#############################################################################
# one executable for each unit test file
#############################################################################

file(GLOB files "src/unit-*.cpp")
list(FILTER files EXCLUDE REGEX "src/unit-modinv.cpp$")

foreach(file ${files})
    get_filename_component(file_basename ${file} NAME_WE)
    string(REGEX REPLACE "unit-([^$]+)" "test-\\1" testcase ${file_basename})

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(_test_opts -Wno-deprecated -Wno-float-equal)
    elseif(CMAKE_COMPILER_IS_GNUCXX)
        set(_test_opts -Wno-deprecated -Wno-float-equal -fconstexpr-depth=30)
    else()
        set(_test_opts "")
    endif()

    ctbignum_add_executable(${testcase}
        SOURCE   ${file}
        CATCH_MAIN
        INCLUDES thirdparty/catch
        DEFINES  CATCH_CONFIG_FAST_COMPILE
        OPTIONS  ${_test_opts})

    # prime_modulus runs Miller-Rabin at compile time for 256-bit moduli
    if(testcase MATCHES "^test-roots")
        lam_target_constexpr_limit(${testcase} CLANG 5000000 GCC 50000000)
    else()
        lam_target_constexpr_limit(${testcase} CLANG 2300000)
    endif()

    add_test(NAME "${testcase}_default"
      COMMAND ${testcase} ${CATCH_TEST_FILTER}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    # Catch v2 exits 0 when no tests run.
    set_tests_properties("${testcase}_default" PROPERTIES
      LABELS "default"
      FAIL_REGULAR_EXPRESSION "No tests ran")

    add_test(NAME "${testcase}_all"
      COMMAND ${testcase} ${CATCH_TEST_FILTER} "*"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    set_tests_properties("${testcase}_all" PROPERTIES
      LABELS "all"
      FAIL_REGULAR_EXPRESSION "No tests ran")

endforeach()
# Check for NTL
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h)
find_library(NTL_LIBRARY NAMES ntl)
find_library(GMP_LIBRARY NAMES gmp)

if(NTL_INCLUDE_DIR AND NTL_LIBRARY AND GMP_LIBRARY)
    ctbignum_add_executable(unit-modinv
        SOURCE   src/unit-modinv.cpp
        CATCH_MAIN
        LIBS     ${NTL_LIBRARY} ${GMP_LIBRARY}
        INCLUDES thirdparty/catch ${NTL_INCLUDE_DIR})
    lam_target_constexpr_limit(unit-modinv CLANG 2300000)

    add_test(NAME unit-modinv_run COMMAND unit-modinv)
    set_tests_properties(unit-modinv_run PROPERTIES
      FAIL_REGULAR_EXPRESSION "No tests ran")
else()
    message(STATUS "NTL or GMP not found. Skipping modinv tests.")
endif()

ctbignum_add_executable(check_kyber_prime SOURCE src/test_kyber_prime.cpp)
