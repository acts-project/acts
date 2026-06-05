# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Declare the core library test(s).
traccc_add_test(core
   "test_algorithm.cpp"
   "test_module_map.cpp"
   "test_pvalue.cpp"
   "particle.cpp"
   LINK_LIBRARIES GTest::gtest_main traccc_tests_common
   traccc::core traccc::io)

if( TRACCC_USE_ROOT )
   find_package( ROOT COMPONENTS Core RIO Hist REQUIRED )
   target_link_libraries( traccc_test_core
      PRIVATE ROOT::Core )
   target_compile_definitions( traccc_test_core
      PRIVATE TRACCC_HAVE_ROOT )
endif()
