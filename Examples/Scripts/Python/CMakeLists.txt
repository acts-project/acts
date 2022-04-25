file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
install(FILES ${files} DESTINATION python)
