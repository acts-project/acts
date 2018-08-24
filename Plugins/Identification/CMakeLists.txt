add_library(IdentificationPlugin INTERFACE)
target_include_directories(IdentificationPlugin INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(IdentificationPlugin INTERFACE Acts::Core)

install(TARGETS IdentificationPlugin EXPORT ActsIdentificationTargets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
