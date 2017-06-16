add_executable (BoundaryCheckTests BoundaryCheckTests.cpp)
target_link_libraries (BoundaryCheckTests PRIVATE ACTSCore)
add_test (NAME BoundaryCheckTests COMMAND BoundaryCheckTests)

add_executable (BoundaryCheckBenchmark BoundaryCheckBenchmark.cpp)
target_link_libraries (BoundaryCheckBenchmark PRIVATE ACTSCore)

add_executable (SurfaceTests SurfaceTests.cpp)
target_link_libraries (SurfaceTests PRIVATE ACTSCore)
add_test (NAME SurfaceTests COMMAND SurfaceTests)

add_executable (LineSurfaceTests LineSurfaceTests.cpp)
target_link_libraries (LineSurfaceTests PRIVATE ACTSCore)
add_test (NAME LineSurfaceTests COMMAND LineSurfaceTests)

add_executable (PerigeeSurfaceTests PerigeeSurfaceTests.cpp)
target_link_libraries (PerigeeSurfaceTests PRIVATE ACTSCore)
add_test (NAME PerigeeSurfaceTests COMMAND PerigeeSurfaceTests)

add_executable (StrawSurfaceTests StrawSurfaceTests.cpp)
target_link_libraries (StrawSurfaceTests PRIVATE ACTSCore)
add_test (NAME StrawSurfaceTests COMMAND StrawSurfaceTests)

add_executable (CylinderSurfaceTests CylinderSurfaceTests.cpp)
target_link_libraries (CylinderSurfaceTests PRIVATE ACTSCore)
add_test (NAME CylinderSurfaceTests COMMAND CylinderSurfaceTests)

add_executable (ConeSurfaceTests ConeSurfaceTests.cpp)
target_link_libraries (ConeSurfaceTests PRIVATE ACTSCore)
add_test (NAME ConeSurfaceTests COMMAND ConeSurfaceTests)

add_executable (DiscSurfaceTests DiscSurfaceTests.cpp)
target_link_libraries (DiscSurfaceTests PRIVATE ACTSCore)
add_test (NAME DiscSurfaceTests COMMAND DiscSurfaceTests)

add_executable (PlaneSurfaceTests PlaneSurfaceTests.cpp)
target_link_libraries (PlaneSurfaceTests PRIVATE ACTSCore)
add_test (NAME PlaneSurfaceTests COMMAND PlaneSurfaceTests)

add_executable (SurfaceBoundsTests SurfaceBoundsTests.cpp)
target_link_libraries (SurfaceBoundsTests PRIVATE ACTSCore)
add_test (NAME SurfaceBoundsTests COMMAND SurfaceBoundsTests)

add_executable (InfiniteBoundsTests InfiniteBoundsTests.cpp)
target_link_libraries (InfiniteBoundsTests PRIVATE ACTSCore)
add_test (NAME InfiniteBoundsTests COMMAND InfiniteBoundsTests)

add_executable (LineBoundsTests LineBoundsTests.cpp)
target_link_libraries (LineBoundsTests PRIVATE ACTSCore)
add_test (NAME LineBoundsTests COMMAND LineBoundsTests)

add_executable (CylinderBoundsTests CylinderBoundsTests.cpp)
target_link_libraries (CylinderBoundsTests PRIVATE ACTSCore)
add_test (NAME CylinderBoundsTests COMMAND CylinderBoundsTests)

add_executable (ConeBoundsTests ConeBoundsTests.cpp)
target_link_libraries (ConeBoundsTests PRIVATE ACTSCore)
add_test (NAME ConeBoundsTests COMMAND ConeBoundsTests)

add_executable (DiamondBoundsTests DiamondBoundsTests.cpp)
target_link_libraries (DiamondBoundsTests PRIVATE ACTSCore)
add_test (NAME DiamondBoundsTests COMMAND DiamondBoundsTests)

add_executable (EllipseBoundsTests EllipseBoundsTests.cpp)
target_link_libraries (EllipseBoundsTests PRIVATE ACTSCore)
add_test (NAME EllipseBoundsTests COMMAND EllipseBoundsTests)

add_executable (TriangleBoundsTests TriangleBoundsTests.cpp)
target_link_libraries (TriangleBoundsTests PRIVATE ACTSCore)
add_test (NAME TriangleBoundsTests COMMAND TriangleBoundsTests)

add_executable (RectangleBoundsTests RectangleBoundsTests.cpp)
target_link_libraries (RectangleBoundsTests PRIVATE ACTSCore)
add_test (NAME RectangleBoundsTests COMMAND RectangleBoundsTests)

add_executable (TrapezoidBoundsTests TrapezoidBoundsTests.cpp)
target_link_libraries (TrapezoidBoundsTests PRIVATE ACTSCore)
add_test (NAME TrapezoidBoundsTests COMMAND TrapezoidBoundsTests)

add_executable (RadialBoundsTests RadialBoundsTests.cpp)
target_link_libraries (RadialBoundsTests PRIVATE ACTSCore)
add_test (NAME RadialBoundsTests COMMAND RadialBoundsTests)

add_executable (DiscTrapezoidalBoundsTests DiscTrapezoidalBoundsTests.cpp)
target_link_libraries (DiscTrapezoidalBoundsTests PRIVATE ACTSCore)
add_test (NAME DiscTrapezoidalBoundsTests COMMAND DiscTrapezoidalBoundsTests)
