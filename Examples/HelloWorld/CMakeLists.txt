acts_add_library(
    ExamplesHelloWorld
    SHARED
    src/HelloLoggerAlgorithm.cpp
    src/HelloRandomAlgorithm.cpp
    src/HelloWhiteBoardAlgorithm.cpp
)

target_include_directories(
    ActsExamplesHelloWorld
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(ActsExamplesHelloWorld PRIVATE Acts::ExamplesFramework)

acts_compile_headers(ExamplesHelloWorld GLOB "include/**/*.hpp")

add_executable(ActsHelloWorldExample src/HelloWorld.cpp)
target_link_libraries(
    ActsHelloWorldExample
    PRIVATE Acts::ExamplesFramework Acts::ExamplesHelloWorld
)

install(
    TARGETS ActsHelloWorldExample
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
