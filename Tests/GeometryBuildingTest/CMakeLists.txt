################################################################################
# Package: GeometryBuildingTest
################################################################################
gaudi_subdir(GeometryBuildingTest v1r0)

gaudi_depends_on_subdirs(GaudiKernel
			 Core/CoreInterfaces
			 Geometry/GeometryInterfaces
			 Geometry/Detector)

add_definitions(-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")
add_definitions(-DACTS_GAUDI)

gaudi_install_headers(GeometryBuildingTest)

gaudi_add_module(GeometryBuildingTest
		  src/*.cxx
		  INCLUDE_DIRS GaudiKernel CoreInterfaces GeometryInterfaces Detector
		  LINK_LIBRARIES GaudiKernel Detector)
