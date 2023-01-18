message(WARNING "The examples binaries are deprecated and may be removed in the future. Consider using the python bindings.")

# shared code
add_subdirectory(Common)
add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)

# tools
add_subdirectory(Digitization)
add_subdirectory(Fatras)
add_subdirectory_if(Geant4 ACTS_BUILD_EXAMPLES_GEANT4)
add_subdirectory(Generators)
add_subdirectory(Geometry)
add_subdirectory(HelloWorld)
add_subdirectory_if(HepMC3 ACTS_BUILD_EXAMPLES_HEPMC3)
add_subdirectory(MagneticField)
add_subdirectory(MaterialMapping)
add_subdirectory(Misc)
add_subdirectory(Propagation)
add_subdirectory(Reconstruction)
add_subdirectory(Show)
add_subdirectory(Vertexing)
add_subdirectory_if(Alignment ACTS_BUILD_ALIGNMENT)
