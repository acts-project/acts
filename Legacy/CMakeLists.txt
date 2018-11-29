add_library(ActsLegacy INTERFACE)

target_include_directories (ActsLegacy SYSTEM INTERFACE ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
target_include_directories(ActsLegacy INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

install(TARGETS ActsLegacy EXPORT ActsLegacyTargets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

