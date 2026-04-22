set(unittest_extra_libraries ActsCore Acts::Alignment)

add_unittest(Alignment AlignmentTests.cpp)

add_unittest(AlignmentHierarchy HierarchyTests.cpp)
