set(unittest_extra_libraries Acts::Fatras)

add_unittest(FatrasBarcode BarcodeTests.cpp)
add_unittest(FatrasHit HitTests.cpp)
add_unittest(FatrasParticle ParticleTests.cpp)
