# manually trigger documentation build
option (BUILD_DOC "build documentation" OFF)

if (BUILD_DOC)
  # we need doxygen
  find_package (Doxygen QUIET)

  if (DOXYGEN_FOUND)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
    add_custom_target (doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
			 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			 COMMENT "Generating API documentation with Doxygen" VERBATIM)
    install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTS.tag DESTINATION doc)
  # do not fail if doxygen is not available, just inform and continue
  else (DOXYGEN_FOUND)
    message (STATUS "doxygen was not found -> skip generation of documentation")
  endif (DOXYGEN_FOUND)
endif (BUILD_DOC)  
