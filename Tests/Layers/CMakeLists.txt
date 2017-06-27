add_executable (LayerTests LayerTests.cpp)
target_link_libraries (LayerTests PRIVATE ACTSCore)
add_test (NAME LayerTests COMMAND LayerTests)

