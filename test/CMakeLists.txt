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

    add_executable(${testcase} $<TARGET_OBJECTS:catch_main> ${file})
    set_target_properties(${testcase} PROPERTIES
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
    )

    target_compile_definitions(${testcase} PRIVATE CATCH_CONFIG_FAST_COMPILE)
    target_include_directories(${testcase} PRIVATE "thirdparty/catch")
    
    # Link against the ctbignum module library
    target_link_libraries(${testcase} PRIVATE lam::ctbignum)

    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        target_compile_options(${testcase} PRIVATE -Wno-deprecated -Wno-float-equal -fconstexpr-steps=2300000)
    elseif(CMAKE_COMPILER_IS_GNUCXX)
        target_compile_options(${testcase} PRIVATE -Wno-deprecated -Wno-float-equal -fconstexpr-depth=30)
    endif()

    add_test(NAME "${testcase}_default"
      COMMAND ${testcase} ${CATCH_TEST_FILTER}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    set_tests_properties("${testcase}_default" PROPERTIES LABELS "default")

    add_test(NAME "${testcase}_all"
      COMMAND ${testcase} ${CATCH_TEST_FILTER} "*"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    set_tests_properties("${testcase}_all" PROPERTIES LABELS "all")

endforeach()
# Check for NTL
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h)
find_library(NTL_LIBRARY NAMES ntl)
find_library(GMP_LIBRARY NAMES gmp)

if(NTL_INCLUDE_DIR AND NTL_LIBRARY AND GMP_LIBRARY)
    add_executable(test-modinv $<TARGET_OBJECTS:catch_main> src/unit-modinv.cpp)
    target_link_libraries(test-modinv PRIVATE lam::ctbignum ${NTL_LIBRARY} ${GMP_LIBRARY})
    target_include_directories(test-modinv PRIVATE "thirdparty/catch" ${NTL_INCLUDE_DIR})
    
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(test-modinv PRIVATE -fprebuilt-module-path=${CMAKE_BINARY_DIR}/CMakeFiles/lam_ctbignum.dir/include/ctbignum)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(test-modinv PRIVATE -fmodules-ts)
    endif()

    # We need to manually add the dependency on the library's module file for the test executable
    add_dependencies(test-modinv lam_ctbignum)

    # Set the C++ standard for the test target
    set_property(TARGET test-modinv PROPERTY CXX_STANDARD 23)
    set_property(TARGET test-modinv PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET test-modinv PROPERTY CXX_EXTENSIONS OFF)

    # Enable module scanning for the test
    set_property(TARGET test-modinv PROPERTY CXX_SCAN_FOR_MODULES ON)
    set_property(TARGET test-modinv PROPERTY CXX_MODULE_STD ON)

    add_test(NAME test-modinv_run COMMAND test-modinv)
else()
    message(STATUS "NTL or GMP not found. Skipping modinv tests.")
endif()
