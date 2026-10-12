find_package(EDM4HEP 1.0 CONFIG REQUIRED)
find_package(podio 1.7 CONFIG REQUIRED)

acts_add_library(
    ExamplesCalorimeterEDM4hep
    src/EDM4hepCalorimeterInputConverter.cpp
    src/EDM4hepCalorimeterOutputConverter.cpp
    src/EDM4hepCalorimeterClusterOutputConverter.cpp
    src/EDM4hepCalorimeterMetadata.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesCalorimeterEDM4hep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsExamplesCalorimeterEDM4hep
    PUBLIC Acts::ExamplesCalorimeter Acts::ExamplesIoPodio EDM4HEP::edm4hep
)
acts_compile_headers(ExamplesCalorimeterEDM4hep GLOB "include/**/*.hpp")

add_executable(ActsExampleCalorimeterEDM4hep EDM4hepCalorimeterExample.cpp)
target_link_libraries(
    ActsExampleCalorimeterEDM4hep
    PRIVATE Acts::ExamplesCalorimeterEDM4hep EDM4HEP::edm4hepDict
)
install(
    TARGETS ActsExampleCalorimeterEDM4hep
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_subdirectory_if(Jets ACTS_BUILD_EXAMPLES_FASTJET)
