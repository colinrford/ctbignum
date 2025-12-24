cmake_minimum_required(VERSION 3.31.6)

project(ctbignum-benchmarks LANGUAGES CXX)

find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp)
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h)
find_library(NTL_LIBRARIES NAMES ntl)
find_package(Threads REQUIRED)

include(FetchContent)
FetchContent_Declare(
  google-benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.8.3
)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(google-benchmark)

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

    add_executable(${benchmark_name} ${file})

    target_link_libraries(${benchmark_name} PRIVATE ctbignum benchmark::benchmark ${NTL_LIBRARIES} ${GMP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    target_include_directories(${benchmark_name} PRIVATE ${NTL_INCLUDES} ${GMP_INCLUDES})

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${benchmark_name} PRIVATE -fprebuilt-module-path=${CMAKE_BINARY_DIR}/CMakeFiles/ctbignum.dir/include/ctbignum)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(${benchmark_name} PRIVATE -fmodules-ts)
    endif()

    add_dependencies(${benchmark_name} ctbignum)

    set_target_properties(${benchmark_name} PROPERTIES
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        CXX_SCAN_FOR_MODULES ON
        CXX_MODULE_STD ON
    )
endforeach()

