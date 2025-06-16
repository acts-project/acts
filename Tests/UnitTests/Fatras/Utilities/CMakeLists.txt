set(unittest_extra_libraries ActsFatras)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    message(STATUS "No Unit test for GammaDistribution on MacOS X")
else()
    add_unittest(FatrasGammaDistribution GammaDistributionTests.cpp)
endif()
