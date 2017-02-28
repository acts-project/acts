add_executable (BoundaryCheckTests BoundaryCheckTests.cpp)
target_link_libraries (BoundaryCheckTests PRIVATE ACTSCore)
add_test (NAME BoundaryCheckTests COMMAND BoundaryCheckTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST BoundaryCheckTests TARGETS BoundaryCheckTests)

add_executable (SurfaceTests SurfaceTests.cpp)
target_link_libraries (SurfaceTests PRIVATE ACTSCore)
add_test (NAME SurfaceTests COMMAND SurfaceTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST SurfaceTests TARGETS SurfaceTests)

add_executable (LineSurfaceTests LineSurfaceTests.cpp)
target_link_libraries (LineSurfaceTests PRIVATE ACTSCore)
add_test (NAME LineSurfaceTests COMMAND LineSurfaceTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST LineSurfaceTests TARGETS LineSurfaceTests)

add_executable (plane_surface_tests plane_surface_tests.cpp)
target_link_libraries (plane_surface_tests PRIVATE ACTSCore)
add_test (NAME plane_surface_tests COMMAND plane_surface_tests)
acts_add_test_to_cdash_project (PROJECT ACore TEST plane_surface_tests TARGETS plane_surface_tests)

add_executable (SurfaceBoundsTests SurfaceBoundsTests.cpp)
target_link_libraries (SurfaceBoundsTests PRIVATE ACTSCore)
add_test (NAME SurfaceBoundsTests COMMAND SurfaceBoundsTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST SurfaceBoundsTests TARGETS SurfaceBoundsTests)

add_executable (RectangleBoundsTests RectangleBoundsTests.cpp)
target_link_libraries (RectangleBoundsTests PRIVATE ACTSCore)
add_test (NAME RectangleBoundsTests COMMAND RectangleBoundsTests)
acts_add_test_to_cdash_project (PROJECT ACore TEST RectangleBoundsTests TARGETS RectangleBoundsTests)