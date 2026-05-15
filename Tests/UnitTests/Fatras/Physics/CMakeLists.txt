set(unittest_extra_libraries Acts::Fatras)

add_unittest(FatrasBetheBloch BetheBlochTests.cpp)
add_unittest(FatrasBetheHeitler BetheHeitlerTests.cpp)
add_unittest(FatrasScattering ScatteringTests.cpp)
add_unittest(FatrasPhotonConversion PhotonConversionTests.cpp)
