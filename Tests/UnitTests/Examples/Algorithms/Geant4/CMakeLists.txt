set(unittest_extra_libraries ActsExamplesGeant4)

add_unittest(Geant4EventContext EventContextTests.cpp)
add_unittest(Geant4SensitiveSurfaceMapper SensitiveSurfaceMapperTests.cpp)
add_unittest(Geant4RunManager RunManagerTests.cpp)
add_unittest(Geant4PreassignedDecay PreassignedDecayTests.cpp)
target_link_libraries(ActsUnitTestGeant4PreassignedDecay PRIVATE HepMC3::HepMC3)

if(Geant4_multithreaded_FOUND)
    add_unittest(Geant4NativeWorkers NativeWorkerTests.cpp)
    target_sources(
        ActsUnitTestGeant4NativeWorkers
        PRIVATE NativeWorkerTransport.cpp NativeWorkerExchangeTests.cpp
    )
    add_test(
        NAME Geant4NativeWorkerReference
        COMMAND
            ActsUnitTestGeant4NativeWorkers
            --run_test=HandOffEventsToNativeWorkers -- reference
    )
    add_test(
        NAME Geant4NativeWorkerShortInput
        COMMAND
            ActsUnitTestGeant4NativeWorkers
            --run_test=HandOffEventsToNativeWorkers -- short
    )
    set_tests_properties(
        Geant4NativeWorkers
        Geant4NativeWorkerShortInput
        PROPERTIES ENVIRONMENT "G4FORCENUMBEROFTHREADS=2" TIMEOUT 60
    )
    set_tests_properties(
        Geant4NativeWorkerReference
        PROPERTIES ENVIRONMENT "G4FORCENUMBEROFTHREADS=1" TIMEOUT 60
    )
endif()

set(unittest_extra_libraries ActsExamplesGeant4 ActsExamplesDetectorTelescope)
add_unittest(Geant4SimulationActions Geant4SimulationActionsTests.cpp)

add_unittest(Geant4Session Geant4SessionTests.cpp)
