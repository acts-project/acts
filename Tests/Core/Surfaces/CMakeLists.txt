add_executable (BoundaryCheckTests BoundaryCheckTests.cpp)
target_link_libraries (BoundaryCheckTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME BoundaryCheckUnitTest COMMAND BoundaryCheckTests)

add_executable (BoundaryCheckBenchmark BoundaryCheckBenchmark.cpp)
target_link_libraries (BoundaryCheckBenchmark PRIVATE ActsCore ActsTestsCommonHelpers)
install (TARGETS BoundaryCheckBenchmark
	     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/Profiling)

add_executable (SurfaceTests SurfaceTests.cpp)
target_link_libraries (SurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME SurfaceUnitTest COMMAND SurfaceTests)

add_executable (LineSurfaceTests LineSurfaceTests.cpp)
target_link_libraries (LineSurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME LineSurfaceUnitTest COMMAND LineSurfaceTests)

add_executable (PerigeeSurfaceTests PerigeeSurfaceTests.cpp)
target_link_libraries (PerigeeSurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME PerigeeSurfaceUnitTest COMMAND PerigeeSurfaceTests)

add_executable (StrawSurfaceTests StrawSurfaceTests.cpp)
target_link_libraries (StrawSurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME StrawSurfaceUnitTest COMMAND StrawSurfaceTests)

add_executable (CylinderSurfaceTests CylinderSurfaceTests.cpp)
target_link_libraries (CylinderSurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME CylinderSurfaceUnitTest COMMAND CylinderSurfaceTests)

add_executable (ConeSurfaceTests ConeSurfaceTests.cpp)
target_link_libraries (ConeSurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME ConeSurfaceUnitTest COMMAND ConeSurfaceTests)

add_executable (DiscSurfaceTests DiscSurfaceTests.cpp)
target_link_libraries (DiscSurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME DiscSurfaceUnitTest COMMAND DiscSurfaceTests)

add_executable (PlaneSurfaceTests PlaneSurfaceTests.cpp)
target_link_libraries (PlaneSurfaceTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME PlaneSurfaceUnitTest COMMAND PlaneSurfaceTests)

add_executable (SurfaceBoundsTests SurfaceBoundsTests.cpp)
target_link_libraries (SurfaceBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME SurfaceBoundsUnitTest COMMAND SurfaceBoundsTests)

add_executable (InfiniteBoundsTests InfiniteBoundsTests.cpp)
target_link_libraries (InfiniteBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME InfiniteBoundsUnitTest COMMAND InfiniteBoundsTests)

add_executable (LineBoundsTests LineBoundsTests.cpp)
target_link_libraries (LineBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME LineBoundsUnitTest COMMAND LineBoundsTests)

add_executable (CylinderBoundsTests CylinderBoundsTests.cpp)
target_link_libraries (CylinderBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME CylinderBoundsUnitTest COMMAND CylinderBoundsTests)

add_executable (ConeBoundsTests ConeBoundsTests.cpp)
target_link_libraries (ConeBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME ConeBoundsUnitTest COMMAND ConeBoundsTests)

add_executable (DiamondBoundsTests DiamondBoundsTests.cpp)
target_link_libraries (DiamondBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME DiamondBoundsUnitTest COMMAND DiamondBoundsTests)

add_executable (EllipseBoundsTests EllipseBoundsTests.cpp)
target_link_libraries (EllipseBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME EllipseBoundsTests COMMAND EllipseBoundsTests)

add_executable (TriangleBoundsTests TriangleBoundsTests.cpp)
target_link_libraries (TriangleBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME TriangleBoundsTests COMMAND TriangleBoundsTests)

add_executable (RectangleBoundsTests RectangleBoundsTests.cpp)
target_link_libraries (RectangleBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME RectangleBoundsTests COMMAND RectangleBoundsTests)

add_executable (TrapezoidBoundsTests TrapezoidBoundsTests.cpp)
target_link_libraries (TrapezoidBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME TrapezoidBoundsTests COMMAND TrapezoidBoundsTests)

add_executable (RadialBoundsTests RadialBoundsTests.cpp)
target_link_libraries (RadialBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME RadialBoundsUnitTest COMMAND RadialBoundsTests)

add_executable (DiscTrapezoidalBoundsTests DiscTrapezoidalBoundsTests.cpp)
target_link_libraries (DiscTrapezoidalBoundsTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME DiscTrapezoidalBoundsUnitTest COMMAND DiscTrapezoidalBoundsTests)

add_executable (SurfaceArrayTests SurfaceArrayTests.cpp)
target_link_libraries (SurfaceArrayTests PRIVATE ActsCore ActsTestsCommonHelpers)
add_test (NAME SurfaceArrayUnitTest COMMAND SurfaceArrayTests)
