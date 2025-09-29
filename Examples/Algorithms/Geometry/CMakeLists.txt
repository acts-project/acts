acts_add_library(ExamplesGeometry SHARED src/VolumeAssociationTest.cpp)

target_include_directories(
    ActsExamplesGeometry
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(ActsExamplesGeometry PUBLIC Acts::ExamplesFramework)
