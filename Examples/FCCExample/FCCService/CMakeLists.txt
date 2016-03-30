################################################################################
# Package: FCCService
################################################################################
gaudi_subdir(FCCService v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Core/CoreInterfaces
                         Plugins/DD4hepPlugins/DD4hepGeometry/DD4hepGeometryInterfaces)

find_package(DD4hep)

gaudi_install_headers(FCCService)

gaudi_add_module(FCCService
                 src/*.cxx
                 INCLUDE_DIRS GaudiKernel DD4hep DD4hepGeometryInterfaces CoreInterfaces
                 LINK_LIBRARIES GaudiKernel DD4hep)
