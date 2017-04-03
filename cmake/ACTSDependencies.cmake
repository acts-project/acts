# let the user overwrite verbosity
set (ACTS_FIND_QUIET OFF CACHE BOOL "find external packages quietly")
if (${ACTS_FIND_QUIET})
  set (ACTS_FIND_QUIET_TEXT "QUIET")
else ()
  set (ACTS_FIND_QUIET_TEXT "")
endif()

# specify version numbers for dependencies of ACTS
set (ACTS_BOOST_VERSION "1.62")
set (ACTS_DOXYGEN_VERSION "1.8.11")
set (ACTS_ROOT_VERSION "6.06.04")
