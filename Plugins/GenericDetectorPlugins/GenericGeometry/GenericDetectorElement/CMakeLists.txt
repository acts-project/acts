################################################################################
# Package: GenericDetectorElement
################################################################################
gaudi_subdir(GenericDetectorElement v1r0)


gaudi_depends_on_subdirs(GaudiKernel
			 Core/CoreIdentifier
			 Core/Algebra
			 Core/CoreUtils
			 Geometry/DetectorElementBase
			 Geometry/Surfaces)

add_definitions(-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")
add_definitions(-DACTS_GAUDI)

find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_add_library(GenericDetectorElement
                  src/*.cxx
                  PUBLIC_HEADERS GenericDetectorElement
                  INCLUDE_DIRS GaudiKernel ${EIGEN3_HOME} CoreIdentifier Algebra CoreUtils DetectorElementBase Surfaces
                  LINK_LIBRARIES Surfaces)
