################################################################################
# Package: GenericGeometryTools
################################################################################
gaudi_subdir(GenericGeometryTools v1r0)

gaudi_depends_on_subdirs(GaudiKernel
			 Core/Algebra
			 Core/CoreInterfaces
			 Geometry/GeometryInterfaces
			 Geometry/Volumes
			 Geometry/Detector
			 Geometry/Material
			 Geometry/GeometryUtils
			 Plugins/GenericDetectorPlugins/GenericGeometry/GenericDetectorElement)

add_definitions(-DATS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")

find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_install_headers(GenericGeometryTools)

gaudi_add_module(GenericGeometryTools
		  src/*.cxx
		  INCLUDE_DIRS GaudiKernel ${EIGEN3_HOME} Algebra CoreInterfaces  GeometryInterfaces Volumes Detector Material GeometryUtils GenericDetectorElement
		  LINK_LIBRARIES GaudiKernel  Volumes Detector Surfaces Material GeometryUtils GenericDetectorElement)
