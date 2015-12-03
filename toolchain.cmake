# Special toolchain file that inherits the same heptools version as the
# used projects.
find_file(inherit_heptools_module InheritHEPTools.cmake)

# FIXME: generator versions must be moved to another file
set(generators_versions
    alpgen      2.1.4
    herwig++    2.7.1      # 2.7.0a not available in LCG_70root6
#    hijing      1.383bs.2
#    lhapdf      6.1.4
#    photos++    3.56
#    powheg-box  r2092
#    pythia6     428.2
    pythia8     186
#    rivet       1.9.0
#    tauola++    1.1.4
#    thepeg      1.9.2p1    # 1.9.0a not available in LCG_70root6
    )

# Process the lines of LCG_generators_*.txt file to extract the
# needed generators (variable generators_versions)
macro(lcg_set_generator name hash version dir)
  #message(STATUS "Processing ${name} ${version} (${dir})")
  list(FIND generators_versions ${name} _gen_idx)
  if(NOT _gen_idx EQUAL -1)
    math(EXPR _gen_vers_idx "${_gen_idx} + 1")
    list(GET generators_versions ${_gen_vers_idx} _gen_vers)
    #message(STATUS "  required ${_gen_vers}")
    if("${version}" STREQUAL "${_gen_vers}")
      set(${name}_config_version ${version} CACHE STRING "Version of ${name}")
      mark_as_advanced(${name}_config_version)
      set(${name}_native_version ${${name}_config_version})
      set(${name}_home ${LCG_releases}/${dir})
      #message(STATUS "  adding ${LCG_releases}/${dir} to search path")
      list(INSERT CMAKE_PREFIX_PATH 0 ${LCG_releases}/${dir})
      list(APPEND found_generators ${name})
    #else()
    #  message(STATUS "  wrong version (${name} ${version})")
    endif()
  #else()
  #  message(STATUS "  ${name} not required")
  endif()
endmacro()

# Check that all the generators in generators_versions are found in the
# LCG_generators_*.txt file.
macro(check_generators)
  set(_idx 0)
  set(_missing_generators)
  list(LENGTH generators_versions _length)
  while(_idx LESS _length)
    list(GET generators_versions ${_idx} _gen_name)
    list(FIND found_generators ${_gen_name} _found_idx)
    if(_found_idx EQUAL -1)
      set(_missing_generators "${_missing_generators} ${_gen_name}")
    endif()
    math(EXPR _idx "${_idx} + 2")
  endwhile()
  if(_missing_generators)
    message(FATAL_ERROR "Missing generators: ${_missing_generators}")
  endif()
endmacro()

# this check is needed because the toolchain seem to be called a second time
# without the proper cache
if(inherit_heptools_module)
  include(${inherit_heptools_module})
  inherit_heptools()

  if(LCG_TOOLCHAIN_INFO)
    string(REPLACE "LCG_externals_" "LCG_generators_" LCG_GENERATOR_INFO "${LCG_TOOLCHAIN_INFO}")
    if(LCG_GENERATOR_INFO STREQUAL LCG_TOOLCHAIN_INFO OR NOT EXISTS ${LCG_GENERATOR_INFO})
      message(FATAL_ERROR "No generators info for ${heptools_version}")
    endif()
    message(STATUS "Using generators infos from ${LCG_GENERATOR_INFO}")
    set(found_generators)
    file(STRINGS ${LCG_GENERATOR_INFO} _lcg_infos)
    foreach(_l ${_lcg_infos})
      if(NOT _l MATCHES "^(PLATFORM|VERSION):")
        string(REGEX REPLACE "; *" ";" _l "${_l}")
        lcg_set_generator(${_l})
      endif()
    endforeach()
    check_generators()
  else()
    # This check is to handle the special case of test builds internal to CMake
    if(NOT CMAKE_SOURCE_DIR MATCHES "CMakeTmp")
      message(FATAL_ERROR "Only LCG >= 68 is supported")
    endif()
  endif()
endif()
