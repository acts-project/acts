################################################################################
# Package: ExtrapolationUtils
################################################################################
gaudi_subdir(ExtrapolationUtils v1r0)


gaudi_depends_on_subdirs(GaudiKernel
			 Core/Algebra
			 EventData/EventDataUtils
			 EventData/TrackParameters
			 EventData/NeutralParameters
			 Geometry/Surfaces
			 Geometry/Material
			 Geometry/Detector
			 Geometry/GeometryUtils
			 MagneticField/MagneticFieldUtils)

add_definitions(-DATS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")

find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_add_library(ExtrapolationUtils
		  src/*.cxx
		  PUBLIC_HEADERS ExtrapolationUtils
		  INCLUDE_DIRS ${EIGEN3_HOME} GaudiKernel Algebra EventDataUtils TrackParameters NeutralParameters Surfaces Material Detector GeometryUtils MagneticFieldUtils
		  LINK_LIBRARIES GaudiKernel TrackParameters NeutralParameters Surfaces Material Detector GeometryUtils MagneticFieldUtils)
