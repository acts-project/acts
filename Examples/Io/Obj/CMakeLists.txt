acts_add_library(
    ExamplesIoObj
    src/ObjTrackingGeometryWriter.cpp
    src/ObjSimHitWriter.cpp
    src/ObjPropagationStepsWriter.cpp
)

target_include_directories(
    ActsExamplesIoObj
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoObj
    PUBLIC Acts::ExamplesFramework Threads::Threads
)

acts_compile_headers(ExamplesIoObj GLOB "include/**/*.hpp")
