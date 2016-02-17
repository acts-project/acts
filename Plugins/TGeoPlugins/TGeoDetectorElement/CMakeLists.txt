################################################################################
# Package: TGeoDetectorElement
################################################################################
gaudi_subdir(TGeoDetectorElement v1r0)

gaudi_depends_on_subdirs(Core/Algebra
			 Core/CoreIdentifier
			 Geometry/DetectorElementBase
			 Geometry/Surfaces)

find_package(Eigen3 REQUIRED)
SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )
find_package(ROOT COMPONENTS Geom)

gaudi_add_library(TGeoDetectorElement
                  src/*.cxx
                  PUBLIC_HEADERS TGeoDetectorElement
                  INCLUDE_DIRS ${EIGEN3_HOME} ROOT DetectorElementBase CoreIdentifier Surfaces
                  LINK_LIBRARIES ROOT Surfaces)
