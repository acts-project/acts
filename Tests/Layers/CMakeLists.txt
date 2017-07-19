add_executable (LayerTests LayerTests.cpp)
target_link_libraries (LayerTests PRIVATE ACTSCore)
add_test (NAME LayerTests COMMAND LayerTests)

add_executable (PlaneLayerTests PlaneLayerTests.cpp)
target_link_libraries (PlaneLayerTests PRIVATE ACTSCore)
add_test (NAME PlaneLayerTests COMMAND PlaneLayerTests)

add_executable (DiscLayerTests DiscLayerTests.cpp)
target_link_libraries (DiscLayerTests PRIVATE ACTSCore)
add_test (NAME DiscLayerTests COMMAND DiscLayerTests)

add_executable (CylinderLayerTests CylinderLayerTests.cpp)
target_link_libraries (CylinderLayerTests PRIVATE ACTSCore)
add_test (NAME CylinderLayerTests COMMAND CylinderLayerTests)

add_executable (ConeLayerTests ConeLayerTests.cpp)
target_link_libraries (ConeLayerTests PRIVATE ACTSCore)
add_test (NAME ConeLayerTests COMMAND ConeLayerTests)


