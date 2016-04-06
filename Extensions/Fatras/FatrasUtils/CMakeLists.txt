################################################################################
# Package: FatrasUtils
################################################################################
gaudi_subdir(FatrasUtils v1r0)

gaudi_depends_on_subdirs(GaudiKernel)

add_definitions(-DACTS_GAUDI)

gaudi_add_library(FatrasUtils
                  src/*.cxx
                  PUBLIC_HEADERS GFatrasUtils
                  INCLUDE_DIRS GaudiKernel
                  LINK_LIBRARIES GaudiKernel)
