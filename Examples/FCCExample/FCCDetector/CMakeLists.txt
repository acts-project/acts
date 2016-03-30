################################################################################
# Package: FCCDetector
################################################################################
gaudi_subdir(FCCDetector v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Plugins/DD4hepPlugins/DD4hepGeometry/DD4hepDetectorElement)

find_package(DD4hep)
find_package(Geant4)
include(${Geant4_USE_FILE})

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )
include( DD4hep )

find_package(ROOT COMPONENTS MathCore GenVector Geom REQUIRED)

gaudi_add_module(FCCDetector
                 src/*.cxx
                 INCLUDE_DIRS  GaudiKernel DD4hep ROOT DetExtensions Geant4 DD4hepDetectorElement
                 LINK_LIBRARIES GaudiKernel DD4hep ROOT Geant4 DD4hepDetectorElement)

set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
message(STATUS "LIBRARY_OUTPUT_PATH -> ${LIBRARY_OUTPUT_PATH}")
dd4hep_generate_rootmap(FCCDetector)

