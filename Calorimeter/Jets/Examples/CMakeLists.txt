acts_add_library(
    ExamplesCalorimeterJets
    src/CalorimeterJetAlgorithm.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)
target_include_directories(
    ActsExamplesCalorimeterJets
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsExamplesCalorimeterJets
    PUBLIC Acts::CalorimeterJets Acts::ExamplesFramework
)
acts_compile_headers(ExamplesCalorimeterJets GLOB "include/**/*.hpp")

add_executable(ActsExampleCalorimeterJets ../../Examples/CalorimeterExample.cpp)
target_compile_definitions(
    ActsExampleCalorimeterJets
    PRIVATE ACTS_EXAMPLE_CALORIMETER_JETS
)
target_link_libraries(
    ActsExampleCalorimeterJets
    PRIVATE Acts::ExamplesCalorimeter Acts::ExamplesCalorimeterJets
)
install(
    TARGETS ActsExampleCalorimeterJets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
