# Special toolchain file that inherits the same heptools version as the
# used projects.
find_file(inherit_heptools_module InheritHEPTools.cmake)

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
#        lcg_set_generator(${_l})
      endif()
    endforeach()
  else()
    # This check is to handle the special case of test builds internal to CMake
    if(NOT CMAKE_SOURCE_DIR MATCHES "CMakeTmp")
      message(FATAL_ERROR "Only LCG >= 68 is supported")
    endif()
  endif()
endif()
