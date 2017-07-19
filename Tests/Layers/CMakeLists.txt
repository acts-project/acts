add_executable (LayerTests LayerTests.cpp)
target_link_libraries (LayerTests PRIVATE ACTSCore)
add_test (NAME LayerTests COMMAND LayerTests)

add_executable (PlaneLayerTests PlaneLayerTests.cpp)
target_link_libraries (PlaneLayerTests PRIVATE ACTSCore)
add_test (NAME PlaneLayerTests COMMAND PlaneLayerTests)

