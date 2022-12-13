add_executable(
  ActsExampleHelloWorld
  HelloWorld.cpp
  HelloLoggerAlgorithm.cpp
  HelloRandomAlgorithm.cpp
  HelloService.cpp
  HelloWhiteBoardAlgorithm.cpp)
target_link_libraries(
  ActsExampleHelloWorld
  PRIVATE
    ActsCore
    ActsExamplesFramework ActsExamplesCommon)

install(
  TARGETS ActsExampleHelloWorld
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
