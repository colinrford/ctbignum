include_guard(GLOBAL)

# ctbignum_add_executable(<name> [SOURCE <file>] [CATCH_MAIN]
#   [LIBS ...] [INCLUDES ...] [DEFINES ...] [OPTIONS ...])
# Executable linked to lam::ctbignum, for test/ and benchmarks/.
function(ctbignum_add_executable NAME)
  cmake_parse_arguments(ARG "CATCH_MAIN" "SOURCE" "LIBS;INCLUDES;DEFINES;OPTIONS" ${ARGN})

  if(NOT ARG_SOURCE)
    set(ARG_SOURCE "${NAME}.cpp")
  endif()

  set(_sources ${ARG_SOURCE})
  if(ARG_CATCH_MAIN)
    list(PREPEND _sources $<TARGET_OBJECTS:catch_main>)
  endif()

  add_executable(${NAME} ${_sources})
  lam_target_warnings(${NAME})

  target_link_libraries(${NAME} PRIVATE lam::ctbignum ${ARG_LIBS})
  if(ARG_INCLUDES)
    target_include_directories(${NAME} SYSTEM PRIVATE ${ARG_INCLUDES})
  endif()
  if(ARG_DEFINES)
    target_compile_definitions(${NAME} PRIVATE ${ARG_DEFINES})
  endif()
  if(ARG_OPTIONS)
    target_compile_options(${NAME} PRIVATE ${ARG_OPTIONS})
  endif()

  # Point the compiler at the library's prebuilt module interface.
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${NAME} PRIVATE
      -fprebuilt-module-path=${CMAKE_BINARY_DIR}/CMakeFiles/lam_ctbignum.dir/include/ctbignum)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(${NAME} PRIVATE -fmodules-ts)
  endif()

  add_dependencies(${NAME} lam_ctbignum)

  set_target_properties(${NAME} PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    CXX_SCAN_FOR_MODULES ON
    CXX_MODULE_STD ON)
endfunction()
