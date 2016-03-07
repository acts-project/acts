################################################################################
# Package: DD4hepGeometryUtils
################################################################################
gaudi_subdir(DD4hepGeometryUtils v1r0)

gaudi_depends_on_subdirs(Core/Algebra
                         Core/CoreInterfaces
                         Geometry/Detector
			 Geometry/Volumes
                         Geometry/Surfaces
			 Geometry/GeometryUtils
			 Geometry/GeometryInterfaces
 			 Plugins/DD4hepPlugins/DD4hepGeometry/DD4hepDetectorElement)

add_definitions(-std=c++14)

find_package(DD4hep)
find_package(ROOT COMPONENTS Geom)
find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_add_library(DD4hepGeometryUtils
                  src/*.cxx
                  PUBLIC_HEADERS DD4hepGeometryUtils 
                  INCLUDE_DIRS DD4hep ROOT ${EIGEN_HOME} Algebra CoreInterfaces  GeometryInterfacesDetector Volumes Surfaces GeometryUtils DD4hepDetectorElement
                  LINK_LIBRARIES DD4hep ROOT Detector Volumes Surfaces DD4hepDetectorElement)
