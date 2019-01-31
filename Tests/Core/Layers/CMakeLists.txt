add_executable (GenericApproachDescriptorTests GenericApproachDescriptorTests.cpp)
target_link_libraries (GenericApproachDescriptorTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME GenericApproachDescriptorUnitTests COMMAND GenericApproachDescriptorTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST GenericApproachDescriptorUnitTests TARGETS GenericApproachDescriptorTests)

add_executable (LayerTests LayerTests.cpp)
target_link_libraries (LayerTests PRIVATE ActsCore)
add_test (NAME LayerUnitTest COMMAND LayerTests)

add_executable (PlaneLayerTests PlaneLayerTests.cpp)
target_link_libraries (PlaneLayerTests PRIVATE ActsCore)
add_test (NAME PlaneLayerUnitTest COMMAND PlaneLayerTests)

add_executable (DiscLayerTests DiscLayerTests.cpp)
target_link_libraries (DiscLayerTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME DiscLayerUnitTest COMMAND DiscLayerTests)

add_executable (CylinderLayerTests CylinderLayerTests.cpp)
target_link_libraries (CylinderLayerTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME CylinderLayerUnitTest COMMAND CylinderLayerTests)

add_executable (ConeLayerTests ConeLayerTests.cpp)
target_link_libraries (ConeLayerTests PRIVATE ActsCore)
add_test (NAME ConeLayerUnitTest COMMAND ConeLayerTests)

add_executable (NavigationLayerTests NavigationLayerTests.cpp)
target_link_libraries (NavigationLayerTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME NavigationLayerUnitTest COMMAND NavigationLayerTests)

add_executable (ProtoLayerTests ProtoLayerTests.cpp)
target_link_libraries (ProtoLayerTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME ProtoLayerUnitTest COMMAND ProtoLayerTests)
