include_guard(GLOBAL)

# lam_target_constexpr_limit(<target> [CLANG <steps>] [GCC <ops>])
# Clang's -fconstexpr-steps and GCC's -fconstexpr-ops-limit count different
# things, so each compiler takes its own number.
function(lam_target_constexpr_limit target)
  cmake_parse_arguments(ARG "" "CLANG;GCC" "" ${ARGN})
  if(ARG_CLANG)
    target_compile_options(${target} PRIVATE
      "$<$<CXX_COMPILER_ID:Clang,AppleClang>:-fconstexpr-steps=${ARG_CLANG}>")
  endif()
  if(ARG_GCC)
    target_compile_options(${target} PRIVATE
      "$<$<CXX_COMPILER_ID:GNU>:-fconstexpr-ops-limit=${ARG_GCC}>")
  endif()
endfunction()
