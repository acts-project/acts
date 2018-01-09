# let the user overwrite verbosity
if (${ACTS_FIND_QUIET})
  set (ACTS_FIND_QUIET_TEXT "QUIET")
else ()
  set (ACTS_FIND_QUIET_TEXT "")
endif()
