add_library(
    ActsExamplesMagneticField
    SHARED
    src/FieldMapRootIo.cpp
    src/FieldMapTextIo.cpp
    src/ScalableBFieldService.cpp
)
target_include_directories(
    ActsExamplesMagneticField
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesMagneticField
    PUBLIC ActsCore ActsExamplesFramework
    PRIVATE std::filesystem
)

acts_compile_headers(ActsExamplesMagneticField GLOB "include/**/*.hpp")

install(
    TARGETS ActsExamplesMagneticField
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
