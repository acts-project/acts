add_subdirectory_if(Json ACTS_BUILD_PLUGIN_JSON)
add_subdirectory_if(Root ACTS_BUILD_EXAMPLES_ROOT)
add_subdirectory(Csv)
add_subdirectory_if(Podio ACTS_BUILD_EXAMPLES_PODIO)

if(TARGET ActsExamplesIoHepMC3)
    set(unittest_extra_libraries ActsExamplesIoHepMC3)
    add_unittest(HepMC3DecayConversion HepMC3/HepMC3DecayConversionTests.cpp)
endif()
