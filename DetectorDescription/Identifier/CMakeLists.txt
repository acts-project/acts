################################################################################
# Package: Identifier
################################################################################
gaudi_subdir(Identifier v1r0)

gaudi_depends_on_subdirs(GaudiKernel)

gaudi_add_library(Identifier
                  src/*.cxx
                  PUBLIC_HEADERS Identifier
                  INCLUDE_DIRS GaudiKernel
                  LINK_LIBRARIES GaudiKernel)
