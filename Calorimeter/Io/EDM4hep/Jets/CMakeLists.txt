acts_add_library(
    ExamplesCalorimeterJetsEDM4hep
    src/EDM4hepCalorimeterJetOutputConverter.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesCalorimeterJetsEDM4hep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src
)
target_link_libraries(
    ActsExamplesCalorimeterJetsEDM4hep
    PUBLIC Acts::ExamplesCalorimeterEDM4hep Acts::ExamplesCalorimeterJets
)
acts_compile_headers(ExamplesCalorimeterJetsEDM4hep GLOB "include/**/*.hpp")

add_executable(
    ActsExampleCalorimeterJetsEDM4hep
    ../EDM4hepCalorimeterExample.cpp
)
target_compile_definitions(
    ActsExampleCalorimeterJetsEDM4hep
    PRIVATE ACTS_CALORIMETER_EDM4HEP_JETS
)
target_link_libraries(
    ActsExampleCalorimeterJetsEDM4hep
    PRIVATE Acts::ExamplesCalorimeterJetsEDM4hep EDM4HEP::edm4hepDict
)
install(
    TARGETS ActsExampleCalorimeterJetsEDM4hep
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
