add_executable (BoundaryCheckTests BoundaryCheckTests.cpp)
target_link_libraries (BoundaryCheckTests PRIVATE ACTSCore)
add_test (NAME BoundaryCheckTests COMMAND BoundaryCheckTests)

add_executable (SurfaceTests SurfaceTests.cpp)
target_link_libraries (SurfaceTests PRIVATE ACTSCore)
add_test (NAME SurfaceTests COMMAND SurfaceTests)

add_executable (LineSurfaceTests LineSurfaceTests.cpp)
target_link_libraries (LineSurfaceTests PRIVATE ACTSCore)
add_test (NAME LineSurfaceTests COMMAND LineSurfaceTests)

add_executable (plane_surface_tests plane_surface_tests.cpp)
target_link_libraries (plane_surface_tests PRIVATE ACTSCore)
add_test (NAME plane_surface_tests COMMAND plane_surface_tests)

add_executable (SurfaceBoundsTests SurfaceBoundsTests.cpp)
target_link_libraries (SurfaceBoundsTests PRIVATE ACTSCore)
add_test (NAME SurfaceBoundsTests COMMAND SurfaceBoundsTests)

add_executable (RectangleBoundsTests RectangleBoundsTests.cpp)
target_link_libraries (RectangleBoundsTests PRIVATE ACTSCore)
add_test (NAME RectangleBoundsTests COMMAND RectangleBoundsTests)
