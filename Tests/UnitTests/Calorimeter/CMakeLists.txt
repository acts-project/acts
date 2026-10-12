set(unittest_extra_libraries Acts::Calorimeter Acts::ExamplesCalorimeter)
add_unittest(CalorimeterResponse CalorimeterResponseTests.cpp)
add_unittest(CalorimeterClustering CalorimeterClusteringTests.cpp)
add_test(NAME CalorimeterExample COMMAND ActsExampleCalorimeter)
add_dependencies(unit_tests ActsExampleCalorimeter)
if(ACTS_BUILD_EXAMPLES_FASTJET)
    list(APPEND unittest_extra_libraries Acts::ExamplesCalorimeterJets)
    add_unittest(CalorimeterJets CalorimeterJetsTests.cpp)
    add_test(NAME CalorimeterJetExample COMMAND ActsExampleCalorimeterJets)
    add_dependencies(unit_tests ActsExampleCalorimeterJets)
endif()
if(ACTS_BUILD_CALORIMETER_EDM4HEP)
    list(APPEND unittest_extra_libraries Acts::ExamplesCalorimeterEDM4hep)
    add_unittest(CalorimeterEDM4hep CalorimeterEDM4hepTests.cpp)
    add_unittest(CalorimeterEDM4hepReconstruction CalorimeterEDM4hepReconstructionTests.cpp)
    target_include_directories(
        ActsUnitTestCalorimeterEDM4hep
        PRIVATE ${PROJECT_SOURCE_DIR}/Calorimeter/Io/EDM4hep
    )
    target_include_directories(
        ActsUnitTestCalorimeterEDM4hepReconstruction
        PRIVATE ${PROJECT_SOURCE_DIR}/Calorimeter/Io/EDM4hep
    )
    add_test(
        NAME CalorimeterEDM4hepExample
        COMMAND
            ActsExampleCalorimeterEDM4hep
            ${CMAKE_CURRENT_BINARY_DIR}/edm4hep-example
    )
    add_dependencies(unit_tests ActsExampleCalorimeterEDM4hep)
    if(ACTS_BUILD_EXAMPLES_FASTJET)
        list(
            APPEND unittest_extra_libraries
            Acts::ExamplesCalorimeterJetsEDM4hep
        )
        add_unittest(CalorimeterEDM4hepJets CalorimeterEDM4hepJetsTests.cpp)
        target_include_directories(
            ActsUnitTestCalorimeterEDM4hepJets
            PRIVATE ${PROJECT_SOURCE_DIR}/Calorimeter/Io/EDM4hep
        )
        add_test(
            NAME CalorimeterEDM4hepJetExample
            COMMAND
                ActsExampleCalorimeterJetsEDM4hep
                ${CMAKE_CURRENT_BINARY_DIR}/edm4hep-jet-example
        )
        add_dependencies(unit_tests ActsExampleCalorimeterJetsEDM4hep)
    endif()
endif()
