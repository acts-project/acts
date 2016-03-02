################################################################################
# Package: DD4hepGeometryUtils
################################################################################
gaudi_subdir(DD4hepGeometryUtils v1r0)

gaudi_depends_on_subdirs(Core/Algebra
                         Core/CoreInterfaces
                         Geometry/Volumes
                         Geometry/Surfaces)

add_definitions(-std=c++14)

find_package(DD4hep)
find_package(ROOT COMPONENTS Geom)
find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_add_library(DD4hepGeometryUtils
                  src/*.cxx
                  PUBLIC_HEADERS DD4hepGeometryUtils 
                  INCLUDE_DIRS DD4hep ROOT ${EIGEN_HOME} Algebra CoreInterfaces Volumes Surfaces
                  LINK_LIBRARIES DD4hep ROOT Volumes Surfaces)
