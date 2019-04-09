add_executable(
    MaterialCompositionTests MaterialCompositionTests
        .cpp) target_include_directories(MaterialCompositionTests PRIVATE ${
    Boost_INCLUDE_DIRS}) target_link_libraries(MaterialCompositionTests PRIVATE ActsCore
                                                   ActsTestsCommonHelpers ${
                                                       Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

    add_executable(MaterialTests MaterialTests
                       .cpp) target_include_directories(MaterialTests PRIVATE ${
        Boost_INCLUDE_DIRS}) target_link_libraries(MaterialTests PRIVATE ActsCore
                                                       ActsTestsCommonHelpers ${
                                                           Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

        add_executable(
            MaterialPropertiesTests MaterialPropertiesTests
                .cpp) target_include_directories(MaterialPropertiesTests PRIVATE ${
            Boost_INCLUDE_DIRS}) target_link_libraries(MaterialPropertiesTests PRIVATE
                                                           ActsCore ActsTestsCommonHelpers
                                                               ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

            add_executable(
                HomogeneousSurfaceMaterialTests HomogeneousSurfaceMaterialTests
                    .cpp) target_include_directories(HomogeneousSurfaceMaterialTests
                                                         PRIVATE ${
                                                             Boost_INCLUDE_DIRS})
                target_link_libraries(
                    HomogeneousSurfaceMaterialTests PRIVATE ActsCore
                        ActsTestsCommonHelpers ${
                            Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

                    add_executable(
                        BinnedSurfaceMaterialTests BinnedSurfaceMaterialTests
                            .cpp) target_include_directories(BinnedSurfaceMaterialTests
                                                                 PRIVATE ${
                                                                     Boost_INCLUDE_DIRS})
                        target_link_libraries(
                            BinnedSurfaceMaterialTests PRIVATE ActsCore ${
                                Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

                            add_executable(SurfaceMaterialProxyTests
                                               SurfaceMaterialProxyTests.cpp)
                                target_include_directories(
                                    SurfaceMaterialProxyTests PRIVATE ${
                                        Boost_INCLUDE_DIRS})
                                    target_link_libraries(
                                        SurfaceMaterialProxyTests PRIVATE ActsCore
                                                                          ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

                                        add_executable(
                                            InterpolatedMaterialMapTests
                                                InterpolatedMaterialMapTests
                                                    .cpp)
                                            target_include_directories(
                                                InterpolatedMaterialMapTests
                                                    PRIVATE ${
                                                        Boost_INCLUDE_DIRS})
                                                target_link_libraries(
                                                    InterpolatedMaterialMapTests
                                                        PRIVATE ActsCore
                                                                ActsTestsCommonHelpers)

                                                    add_test(
                                                        NAME MaterialUnitTest COMMAND
                                                                              MaterialCompositionTests
                                                                MaterialTests MaterialPropertiesTests
                                                                                HomogeneousSurfaceMaterialTests
                                                                                BinnedSurfaceMaterialTests
                                                                                SurfaceMaterialProxyTests
                                                                                InterpolatedMaterialMapTests)
                                                        acts_add_test_to_cdash_project(
                                                            PROJECT ACore TEST
                                                                MaterialUnitTest TARGETS MaterialCompositionTests
                                                                    MaterialTests MaterialPropertiesTests
                                                                                    HomogeneousSurfaceMaterialTests
                                                                                    BinnedSurfaceMaterialTests
                                                                                    SurfaceMaterialProxyTests
                                                                                    InterpolatedMaterialMapTests)
