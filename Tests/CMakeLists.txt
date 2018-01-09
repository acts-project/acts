# should reuse version defined in main CMakeLists.txt
find_package (Boost REQUIRED COMPONENTS unit_test_framework)
find_package (PythonInterp 2 EXACT REQUIRED)

add_subdirectory (Detector)
add_subdirectory (EventData)
add_subdirectory (Extrapolation)
add_subdirectory (Fitter)
add_subdirectory (Layers)
add_subdirectory (MagneticField)
add_subdirectory (Material)
add_subdirectory_if (MaterialMapping BUILD_MATERIAL_PLUGIN)
add_subdirectory (Propagator)
add_subdirectory (Seeding)
add_subdirectory (Surfaces)
add_subdirectory (Utilities)
add_subdirectory (Volumes)
add_subdirectory (Tools)

# generate coverage report via gcovr
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/coverage)
add_custom_target(
  coverage
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/cmake/gcovr -r . -e '../Tests/' --html --html-details -o coverage/index.html
  COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/cmake/gcovr -r . -e '../Tests/' --branches
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
