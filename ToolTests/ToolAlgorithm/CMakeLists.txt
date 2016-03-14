################################################################################
# Package: ToolAlgorithm
################################################################################
gaudi_subdir(ToolAlgorithm v1r0)

gaudi_depends_on_subdirs(GaudiKernel
			 Core/CoreInterfaces
			 ToolTests/ToolTestInterfaces)

gaudi_install_headers(ToolAlgorithm)

gaudi_add_module(ToolAlgorithm
		  src/*.cxx
		  INCLUDE_DIRS GaudiKernel CoreInterfaces ToolTestInterfaces
		  LINK_LIBRARIES GaudiKernel)
