# helper function for parsing function arguments
include (CMakeParseArguments)

# This function adds targets to a CDash (sub)project.
#
# It takes the name of the CDash (sub)project and a list of targets. A
# custom target for the CDash (sub)project is created with the name
# '<PROJECT>Project' if it does not yet exist. All specified targets are
# added as dependencies to this custom target such that one can build the
# whole (sub)project using one target.
# The specified targets and their associated source files have their LABELS
# property set such that it includes the CDash (sub)project name.
#
# Usage: acts_add_targets_to_cdash_project(PROJECT project TARGETS target1 [target2 ...])
function (acts_add_targets_to_cdash_project)
  # get function arguments
  cmake_parse_arguments (ARG "" "PROJECT" "TARGETS" ${ARGN})
  # get (and possibly created) CDash helper target
  set (SUBPROJECT_TARGET "${ARG_PROJECT}Project")
  if (NOT TARGET ${SUBPROJECT_TARGET})
    add_custom_target (${SUBPROJECT_TARGET})
  endif ()
  
  # process all targets
  foreach (_target ${ARG_TARGETS})
    # include in CDash target and assign label
    add_dependencies (${SUBPROJECT_TARGET} ${_target})
    set_property (TARGET ${_target} APPEND PROPERTY LABELS "${ARG_PROJECT}")
    # set label on all source files (needed for per-subproject test coverage)
    get_property (_sources TARGET ${_target} PROPERTY SOURCES)
    foreach (_source_file ${_sources})
      set_property (SOURCE ${_source_file} APPEND PROPERTY LABELS "${ARG_PROJECT}")
    endforeach ()
  endforeach ()
endfunction ()

# This function adds tests to a CDash (sub)project.
#
# It takes the name of the CDash (sub)project and the name of test. Optionally,
# a list of targets which need to be build before the test can run can be
# provided.
# The test has its LABELS property updated such that it includes the CDash
# (sub)project name. If a list of targets is specfied, it is added as dependencies
# to the custom '<PROJECT>Project' target (which is created if it does not yet
# exist).
#
# Usage: acts_add_test_to_cdash_project(PROJECT project TEST test [TARGETS target1 [target2 ...]])
function (acts_add_test_to_cdash_project)
  cmake_parse_arguments (ARG "" "PROJECT;TEST" "TARGETS" ${ARGN})
  set_property (TEST ${ARG_TEST} APPEND PROPERTY LABELS "${ARG_PROJECT}")
  # add target to CDash project
  foreach (_target ${ARG_TARGETS})
    acts_add_targets_to_cdash_project (PROJECT ${ARG_PROJECT} TARGETS ${_target})
  endforeach ()
endfunction ()

