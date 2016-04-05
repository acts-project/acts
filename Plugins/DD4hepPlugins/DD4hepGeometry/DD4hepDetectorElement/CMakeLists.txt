################################################################################
# Package: DD4hepDetectorElement
################################################################################
gaudi_subdir(DD4hepDetectorElement v1r0)

gaudi_depends_on_subdirs(Plugins/TGeoPlugins/TGeoDetectorElement)

find_package(DD4hep)

gaudi_add_library(DD4hepDetectorElement
                  src/*.cxx
                  PUBLIC_HEADERS DD4hepDetectorElement
                  INCLUDE_DIRS DD4hep TGeoDetectorElement
                  LINK_LIBRARIES DD4hep TGeoDetectorElement)
