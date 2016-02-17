################################################################################
# Package: TGeoPlugins
################################################################################
gaudi_subdir(TGeoPlugins v1r0)

gaudi_depends_on_subdirs(Core/Algebra
			 Core/CoreIdentifier
			 Geometry/DetectorElementBase
			 Geometry/Surfaces)

find_package(Eigen3 REQUIRED)
SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )
find_package(ROOT COMPONENTS Geom)

gaudi_add_library(TGeoPlugins
                  src/*.cxx
                  PUBLIC_HEADERS TGeoPlugins
                  INCLUDE_DIRS ${EIGEN3_HOME} ROOT DetectorElementBase CoreIdentifier Surfaces
                  LINK_LIBRARIES ROOT Surfaces)
