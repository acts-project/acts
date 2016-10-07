# allow to turn of unit tests manually
option (BUILD_TESTS "build unit tests" ON)

if (BUILD_TESTS)
   # we need boost
   find_package (Boost 1.60 REQUIRED)

   if (Boost_FOUND)
      message (STATUS "enable ACTS unit tests")
      # delegate
      add_subdirectory (EventData)
      add_subdirectory (Fitter)
      add_subdirectory (Seeding)
      add_subdirectory (Surfaces)
      add_subdirectory (Utilities)
   else (Boost_FOUND)
      message (STATUS "${Red}Boost not found -> disable ACTS unit tests${ColorReset}")
   endif (Boost_FOUND)
else (BUILD_TESTS)
  message (STATUS "${Blue}disable ACTS unit tests${ColorReset}")
endif (BUILD_TESTS)
