find_package (Doxygen)

if (DOXYGEN_FOUND)
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
  add_custom_target (doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
			 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			 COMMENT "Generating API documentation with Doxygen" VERBATIM)
  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc)
  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/ACTS.tag DESTINATION doc)
endif (DOXYGEN_FOUND)
