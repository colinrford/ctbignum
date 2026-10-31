cmake_minimum_required(VERSION 3.31.6)

project(ctbignum-benchmarks LANGUAGES CXX)

## standalone-*: no external dependencies
file(GLOB standalone_files "src/standalone-*.cpp")
foreach(file ${standalone_files})
    get_filename_component(file_basename ${file} NAME_WE)
    ctbignum_add_executable(${file_basename} SOURCE ${file} LIBS ${CMAKE_THREAD_LIBS_INIT})
    if(file_basename MATCHES "roots")
        lam_target_constexpr_limit(${file_basename} CLANG 5000000 GCC 50000000)
    endif()
endforeach()

## bench-*, bench.cpp: upstream comparisons against NTL; skipped if deps are missing
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp)
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h)
find_library(NTL_LIBRARIES NAMES ntl)
find_package(benchmark CONFIG QUIET)

# The config package doesn't say which standard library benchmark was built
# against; a libc++ build fails to link under GCC/libstdc++, so test a link.
if(benchmark_FOUND)
    include(CheckCXXSourceCompiles)
    block()
        # try_compile inherits import std and fails without it being relevant here
        set(CMAKE_CXX_MODULE_STD OFF)
        set(CMAKE_REQUIRED_LIBRARIES benchmark::benchmark)
        set(CMAKE_REQUIRED_QUIET ON)
        check_cxx_source_compiles([=[
            #include <benchmark/benchmark.h>
            int main(int argc, char** argv) {
                benchmark::Initialize(&argc, argv);
                benchmark::RunSpecifiedBenchmarks();
            }
        ]=] CTBIGNUM_BENCHMARK_LINKS)
    endblock()
    if(NOT CTBIGNUM_BENCHMARK_LINKS)
        message(WARNING "google-benchmark at ${benchmark_DIR} doesn't link with "
            "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (likely built "
            "against a different standard library).")
    endif()
endif()

set(_missing "")
if(NOT (GMP_INCLUDE_DIR AND GMP_LIBRARIES))
    list(APPEND _missing GMP)
endif()
if(NOT (NTL_INCLUDE_DIR AND NTL_LIBRARIES))
    list(APPEND _missing NTL)
endif()
if(NOT benchmark_FOUND)
    list(APPEND _missing google-benchmark)
elseif(NOT CTBIGNUM_BENCHMARK_LINKS)
    list(APPEND _missing "google-benchmark (doesn't link)")
endif()

if(NOT _missing)
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
        if(file_basename MATCHES "roots")
            lam_target_constexpr_limit(${benchmark_name} CLANG 5000000 GCC 50000000)
        endif()
    endforeach()
else()
    list(JOIN _missing ", " _missing)
    message(STATUS "Skipping bench-* benchmarks: ${_missing}")
endif()
