add_library(
    ActsExamplesIoHepMC3
    SHARED
    src/HepMC3Reader.cpp
    src/HepMC3Writer.cpp
    src/HepMC3InputConverter.cpp
    src/HepMC3OutputConverter.cpp
    src/HepMC3Util.cpp
)
target_include_directories(
    ActsExamplesIoHepMC3
    SYSTEM
    PUBLIC ${HEPMC3_INCLUDE_DIR}
)
target_include_directories(
    ActsExamplesIoHepMC3
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoHepMC3
    PUBLIC ActsCore ActsExamplesFramework ${HEPMC3_LIBRARIES}
    PRIVATE ActsFatras
)

# The code below conditionally tries to find a number of compression libraries,
# and will then configure the HepMC3 IO target to make use of them
# conditionally.
target_compile_definitions(ActsExamplesIoHepMC3 PRIVATE HEPMC3_USE_COMPRESSION)

message(CHECK_START "Examples HepMC3 IO: checking for zlib")
find_package(ZLIB)
if(ZLIB_FOUND)
    message(CHECK_PASS "found")
    target_link_libraries(ActsExamplesIoHepMC3 PRIVATE ZLIB::ZLIB)
    target_compile_definitions(ActsExamplesIoHepMC3 PRIVATE HEPMC3_Z_SUPPORT=1)
else()
    message(CHECK_PASS "not found")
endif()

message(CHECK_START "Examples HepMC3 IO: checking for lzma")
find_package(LibLZMA)
if(LIBLZMA_FOUND)
    message(CHECK_PASS "found")
    target_link_libraries(ActsExamplesIoHepMC3 PRIVATE LibLZMA::LibLZMA)
    target_compile_definitions(
        ActsExamplesIoHepMC3
        PRIVATE HEPMC3_LZMA_SUPPORT=1
    )
else()
    message(CHECK_PASS "not found")
endif()

message(CHECK_START "Examples HepMC3 IO: checking for bzip2")
find_package(BZip2)
if(BZIP2_FOUND)
    message(CHECK_PASS "found")
    target_link_libraries(ActsExamplesIoHepMC3 PRIVATE BZip2::BZip2)
    target_compile_definitions(
        ActsExamplesIoHepMC3
        PRIVATE HEPMC3_BZ2_SUPPORT=1
    )
else()
    message(CHECK_PASS "not found")
endif()

message(CHECK_START "Examples HepMC3 IO: checking for zstd")
find_package(zstd)
if(zstd_FOUND)
    message(CHECK_PASS "found")
    target_link_libraries(ActsExamplesIoHepMC3 PRIVATE zstd::libzstd)
    target_compile_definitions(
        ActsExamplesIoHepMC3
        PRIVATE HEPMC3_ZSTD_SUPPORT=1
    )
else()
    message(CHECK_PASS "not found")
endif()

install(
    TARGETS ActsExamplesIoHepMC3
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
