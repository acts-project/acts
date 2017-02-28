# allow to turn of unit tests manually
option (BUILD_TESTS "build unit tests" ON)

if (BUILD_TESTS)
   # we need boost
   find_package (Boost ${ACTS_BOOST_VERSION} REQUIRED)

   if (Boost_FOUND)
      message (STATUS "enable ACTS unit tests")
      # delegate
      add_subdirectory (Detector)
      add_subdirectory (EventData)
      add_subdirectory (Fitter)
      add_subdirectory (MagneticField)
      add_subdirectory (Material)
      add_subdirectory (Propagator)
      add_subdirectory (Seeding)
      add_subdirectory (Surfaces)
      add_subdirectory (Utilities)
      add_subdirectory (Volumes)
      add_subdirectory (profiling)
      add_subdirectory (Tools)
   else (Boost_FOUND)
      message (STATUS "${Red}Boost not found -> disable ACTS unit tests${ColorReset}")
   endif (Boost_FOUND)
else (BUILD_TESTS)
  message (STATUS "${Blue}disable ACTS unit tests${ColorReset}")
endif (BUILD_TESTS)
