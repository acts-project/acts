# shared code
add_subdirectory(Common)

# tools
add_subdirectory(EventGenerator)
add_subdirectory_if(Fatras ACTS_BUILD_EXAMPLES_PYTHIA8)
add_subdirectory(Geometry)
add_subdirectory(HelloWorld)
add_subdirectory_if(HepMC3 ACTS_BUILD_EXAMPLES_HEPMC3)
add_subdirectory(MagneticField)
add_subdirectory(MaterialMapping)
add_subdirectory(Misc)
add_subdirectory(Propagation)
add_subdirectory(ReadCsv)
add_subdirectory(Reconstruction)
add_subdirectory_if(Vertexing ACTS_BUILD_EXAMPLES_PYTHIA8)
