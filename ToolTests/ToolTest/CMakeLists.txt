################################################################################
# Package: ToolTest
################################################################################
gaudi_subdir(ToolTest v1r0)

gaudi_depends_on_subdirs(GaudiKernel
			 Core/CoreInterfaces
			 ToolTests/ToolTestInterfaces)

gaudi_install_headers(ToolTest)

gaudi_add_module(ToolTest
                 src/*.cxx
                 INCLUDE_DIRS GaudiKernel CoreInterfaces ToolTestInterfaces 
                  LINK_LIBRARIES GaudiKernel)

