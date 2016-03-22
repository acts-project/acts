################################################################################
# Package: MagneticFieldServices
################################################################################
gaudi_subdir(MagneticFieldServices v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Core/CoreInterfaces
                         MagneticField/MagneticFieldInterfaces)

gaudi_install_headers(MagneticFieldServices)

gaudi_add_module(MagneticFieldServices
                src/*.cxx
                INCLUDE_DIRS GaudiKernel CoreInterfaces MagneticFieldInterfaces
                LINK_LIBRARIES GaudiKernel)
