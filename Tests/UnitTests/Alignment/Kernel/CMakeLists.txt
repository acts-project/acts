set(unittest_extra_libraries Acts::Alignment)

add_unittest(Alignment AlignmentTests.cpp)

add_unittest(AlignmentHierarchy HierarchyTests.cpp)
