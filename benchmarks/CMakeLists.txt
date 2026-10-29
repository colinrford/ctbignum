cmake_minimum_required(VERSION 3.31.6)

project(ctbignum-benchmarks LANGUAGES CXX)

## standalone-*: no external dependencies
file(GLOB standalone_files "src/standalone-*.cpp")
foreach(file ${standalone_files})
    get_filename_component(file_basename ${file} NAME_WE)
    ctbignum_add_executable(${file_basename} SOURCE ${file} LIBS ${CMAKE_THREAD_LIBS_INIT})
endforeach()

## bench-*, bench.cpp: upstream comparisons against NTL; skipped if deps are missing
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp)
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h)
find_library(NTL_LIBRARIES NAMES ntl)
find_package(benchmark CONFIG QUIET)

if(GMP_INCLUDE_DIR AND GMP_LIBRARIES
   AND NTL_INCLUDE_DIR AND NTL_LIBRARIES
   AND benchmark_FOUND)
    find_package(Threads REQUIRED)

    set(GMP_INCLUDES ${GMP_INCLUDE_DIR})
    set(NTL_INCLUDES ${NTL_INCLUDE_DIR})

    file(GLOB files "src/bench-*.cpp")
    list(APPEND files "src/bench.cpp")

    foreach(file ${files})
        get_filename_component(file_basename ${file} NAME_WE)
        string(REGEX REPLACE "bench-([^$]+)" "benchmark-\\1" benchmark_name ${file_basename})
        if("${file_basename}" STREQUAL "bench")
            set(benchmark_name "cbn-benchmark")
        endif()

        ctbignum_add_executable(${benchmark_name}
            SOURCE   ${file}
            LIBS     benchmark::benchmark ${NTL_LIBRARIES} ${GMP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
            INCLUDES ${NTL_INCLUDES} ${GMP_INCLUDES})
    endforeach()
else()
    message(STATUS "NTL, GMP or google-benchmark not found. Skipping bench-* benchmarks.")
endif()
