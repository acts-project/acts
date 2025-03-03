add_library(ActsExamplesITkHelpers SHARED src/ITkDetectorElement.cpp)

target_link_libraries(ActsExamplesITkHelpers PUBLIC ActsCore)

target_include_directories(
    ActsExamplesITkHelpers
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

install(
    TARGETS ActsExamplesITkHelpers
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
