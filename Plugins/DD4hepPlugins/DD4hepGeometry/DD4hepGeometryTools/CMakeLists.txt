################################################################################
# Package: DD4hepGeometryTools
################################################################################
gaudi_subdir(DD4hepGeometryTools v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Core/CoreInterfaces
                         Core/Algebra
                         Geometry/GeometryInterfaces
                         Geometry/Detector
                         Plugins/DD4hepPlugins/DD4hepGeometry/DD4hepGeometryInterfaces
                         Plugins/DD4hepPlugins/DD4hepGeometry/DD4hepGeometryUtils
                         Plugins/DD4hepPlugins/DD4hepGeometry/DD4hepDetectorElement)

add_definitions(-std=c++14)

find_package(DD4hep)
find_package(ROOT COMPONENTS Geom)
find_package(Eigen3 REQUIRED)
SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_install_headers(DD4hepGeometryTools)

gaudi_add_module(DD4hepGeometryTools
                 src/*.cxx
                 INCLUDE_DIRS GaudiKernel ROOT DD4hep ${EIGEN3_HOME} CoreInterfaces Algebra  GeometryInterfaces Detector DD4hepGeometryInterfaces DD4hepGeometryUtils DD4hepDetectorElement
                 LINK_LIBRARIES GaudiKernel ROOT DD4hep Detector DD4hepGeometryUtils DD4hepDetectorElement)
