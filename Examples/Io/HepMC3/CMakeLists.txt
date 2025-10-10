acts_add_library(
    ExamplesIoHepMC3
    SHARED
    src/HepMC3Reader.cpp
    src/HepMC3Writer.cpp
    src/HepMC3InputConverter.cpp
    src/HepMC3OutputConverter.cpp
    src/HepMC3Util.cpp
    src/HepMC3Metadata.cpp
    src/HepMC3FileChain.cpp
    src/HepMC3Normalizer.cpp
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
    PUBLIC Acts::ExamplesFramework ${HEPMC3_LIBRARIES}
    PRIVATE Acts::Fatras nlohmann_json::nlohmann_json
)

# Standalone HepMC3 normalizer executable
find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(hepmc_normalize HepMC3Normalize.cpp)

target_link_libraries(
    hepmc_normalize
    PRIVATE
        Acts::ExamplesIoHepMC3
        Boost::program_options
        nlohmann_json::nlohmann_json
)

# HepMC3 version < 3.3.0 had non-inline pseudo static members that give linker errors:
# https://gitlab.cern.ch/hepmc/HepMC3/-/merge_requests/284
if(${HEPMC3_VERSION} VERSION_GREATER_EQUAL 3.2.7)
    # The code below conditionally tries to find a number of compression libraries,
    # and will then configure the HepMC3 IO targets to make use of them
    # conditionally.
    target_compile_definitions(
        ActsExamplesIoHepMC3
        PRIVATE HEPMC3_USE_COMPRESSION
    )

    message(CHECK_START "Examples HepMC3 IO: checking for zlib")
    find_package(ZLIB QUIET)
    if(ZLIB_FOUND)
        message(CHECK_PASS "found")
        target_link_libraries(ActsExamplesIoHepMC3 PRIVATE ZLIB::ZLIB)
        target_compile_definitions(
            ActsExamplesIoHepMC3
            PRIVATE HEPMC3_Z_SUPPORT=1
        )
        target_link_libraries(hepmc_normalize PRIVATE ZLIB::ZLIB)
        target_compile_definitions(hepmc_normalize PRIVATE HEPMC3_Z_SUPPORT=1)
    else()
        message(CHECK_PASS "not found")
    endif()

    message(CHECK_START "Examples HepMC3 IO: checking for lzma")
    find_package(LibLZMA QUIET)
    if(LIBLZMA_FOUND)
        message(CHECK_PASS "found")
        target_link_libraries(ActsExamplesIoHepMC3 PRIVATE LibLZMA::LibLZMA)
        target_compile_definitions(
            ActsExamplesIoHepMC3
            PRIVATE HEPMC3_LZMA_SUPPORT=1
        )
        target_link_libraries(hepmc_normalize PRIVATE LibLZMA::LibLZMA)
        target_compile_definitions(
            hepmc_normalize
            PRIVATE HEPMC3_LZMA_SUPPORT=1
        )
    else()
        message(CHECK_PASS "not found")
    endif()

    message(CHECK_START "Examples HepMC3 IO: checking for bzip2")
    find_package(BZip2 QUIET)
    if(BZIP2_FOUND)
        message(CHECK_PASS "found")
        target_link_libraries(ActsExamplesIoHepMC3 PRIVATE BZip2::BZip2)
        target_compile_definitions(
            ActsExamplesIoHepMC3
            PRIVATE HEPMC3_BZ2_SUPPORT=1
        )
        target_link_libraries(hepmc_normalize PRIVATE BZip2::BZip2)
        target_compile_definitions(hepmc_normalize PRIVATE HEPMC3_BZ2_SUPPORT=1)
    else()
        message(CHECK_PASS "not found")
    endif()

    message(CHECK_START "Examples HepMC3 IO: checking for zstd")
    find_package(ZSTD QUIET)
    if(ZSTD_FOUND)
        message(CHECK_PASS "found")
        target_link_libraries(ActsExamplesIoHepMC3 PRIVATE ZSTD::ZSTD)
        target_compile_definitions(
            ActsExamplesIoHepMC3
            PRIVATE HEPMC3_ZSTD_SUPPORT=1
        )
        target_link_libraries(hepmc_normalize PRIVATE ZSTD::ZSTD)
        target_compile_definitions(
            hepmc_normalize
            PRIVATE HEPMC3_ZSTD_SUPPORT=1
        )
    else()
        message(CHECK_PASS "not found")
    endif()

    target_compile_definitions(hepmc_normalize PRIVATE HEPMC3_USE_COMPRESSION)
endif()

message(CHECK_START "Examples HepMC3 IO: checking for ROOT IO")
if(TARGET HepMC3::rootIO)
    message(CHECK_PASS "found")
    target_link_libraries(ActsExamplesIoHepMC3 PRIVATE HepMC3::rootIO)
    target_compile_definitions(
        ActsExamplesIoHepMC3
        PRIVATE ACTS_HEPMC3_ROOT_SUPPORT=1
    )
    target_link_libraries(hepmc_normalize PRIVATE HepMC3::rootIO)
    target_compile_definitions(
        hepmc_normalize
        PRIVATE ACTS_HEPMC3_ROOT_SUPPORT=1
    )
else()
    message(CHECK_PASS "not found")
endif()

install(TARGETS hepmc_normalize RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
