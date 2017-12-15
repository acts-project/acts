if (BUILD_DOC)
  # we need doxygen
  find_package (Doxygen ${ACTS_DOXYGEN_VERSION} ${ACTS_FIND_QUIET_TEXT})

  if (DOXYGEN_FOUND)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/figures/ACTSlogo.gif ${CMAKE_CURRENT_BINARY_DIR}/ACTSlogo.gif COPYONLY)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml COPYONLY)
    add_custom_target (doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
			           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			           COMMENT "Generating API documentation with Doxygen" VERBATIM)
    install (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/figures DESTINATION ${CMAKE_INSTALL_DOCDIR}/html/doc OPTIONAL)
    install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR} OPTIONAL)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTS.tag DESTINATION ${CMAKE_INSTALL_DOCDIR} OPTIONAL)
  # do not fail if doxygen is not available, just inform and continue
  else (DOXYGEN_FOUND)
    message (STATUS "${Red}doxygen was not found -> skip generation of documentation${ColorReset}")
  endif (DOXYGEN_FOUND)
else (BUILD_DOC)
  message (STATUS "${Blue}disable generation of doxygen documentation${ColorReset}")
endif (BUILD_DOC)
