acts_add_library(ExamplesFrameworkML SHARED src/NeuralCalibrator.cpp)

target_include_directories(
    ActsExamplesFrameworkML
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesFrameworkML
    PUBLIC Acts::ExamplesFramework Acts::PluginOnnx
)

install(
    TARGETS ActsExamplesFrameworkML
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/ActsExamples DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
