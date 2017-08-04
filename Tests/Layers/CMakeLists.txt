add_executable (GenericApproachDescriptorTests GenericApproachDescriptorTests.cpp)
target_link_libraries (GenericApproachDescriptorTests PRIVATE ACTSCore)
add_test (NAME GenericApproachDescriptorTests COMMAND GenericApproachDescriptorTests)

add_executable (LayerTests LayerTests.cpp)
target_link_libraries (LayerTests PRIVATE ACTSCore)
add_test (NAME LayerUnitTest COMMAND LayerTests)

add_executable (PlaneLayerTests PlaneLayerTests.cpp)
target_link_libraries (PlaneLayerTests PRIVATE ACTSCore)
add_test (NAME PlaneLayerUnitTest COMMAND PlaneLayerTests)

add_executable (DiscLayerTests DiscLayerTests.cpp)
target_link_libraries (DiscLayerTests PRIVATE ACTSCore)
add_test (NAME DiscLayerUnitTest COMMAND DiscLayerTests)

add_executable (CylinderLayerTests CylinderLayerTests.cpp)
target_link_libraries (CylinderLayerTests PRIVATE ACTSCore)
add_test (NAME CylinderLayerUnitTest COMMAND CylinderLayerTests)

add_executable (ConeLayerTests ConeLayerTests.cpp)
target_link_libraries (ConeLayerTests PRIVATE ACTSCore)
add_test (NAME ConeLayerUnitTest COMMAND ConeLayerTests)

add_executable (NavigationLayerTests NavigationLayerTests.cpp)
target_link_libraries (NavigationLayerTests PRIVATE ACTSCore)
add_test (NAME NavigationLayerUnitTest COMMAND NavigationLayerTests)
