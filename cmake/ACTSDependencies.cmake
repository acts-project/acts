# let the user overwrite verbosity
if (${ACTS_FIND_QUIET})
  set (ACTS_FIND_QUIET_TEXT "QUIET")
else ()
  set (ACTS_FIND_QUIET_TEXT "")
endif()

# specify version numbers for dependencies of ACTS
set (ACTS_BOOST_VERSION "1.62")
set (ACTS_EIGEN_VERSION "3.2.9")
set (ACTS_DOXYGEN_VERSION "1.8.11")
set (ACTS_ROOT_VERSION "6.08.00")
set (ACTS_DD4HEP_VERSION "1.02")
