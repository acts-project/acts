find_package (ROOT ${ACTS_ROOT_VERSION} ${ACTS_FIND_QUIET_TEXT} COMPONENTS RIO Tree MathCore)
if (ROOT_FOUND)

  add_executable (PropagationTests PropagationTests.cpp)
  target_link_libraries (PropagationTests PRIVATE ACTSCore)
  add_integration_test(PropagationTests)

else (ROOT_FOUND)
  message (STATUS "${Red}ROOT was not found -> skip build of Integration Tests${ColorReset}")
endif (ROOT_FOUND)
