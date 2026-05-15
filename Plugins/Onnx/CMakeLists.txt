acts_add_library(
    PluginOnnx
    # source files
    src/OnnxRuntimeBase.cpp
    src/MLTrackClassifier.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginOnnx
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginOnnx PUBLIC Acts::Core onnxruntime::onnxruntime)

acts_compile_headers(PluginOnnx GLOB include/**/*.hpp)
