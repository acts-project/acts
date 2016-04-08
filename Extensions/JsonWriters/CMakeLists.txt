################################################################################
# Package: JsonWriters
################################################################################
gaudi_subdir(JsonWriters v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Core/CoreInterfaces
                         Core/Algebra
                         Geometry/GeometryInterfaces
                         Tests/TestInterfaces
                         EventData/ParametersBase
                         Geometry/Surfaces
                         Geometry/Volumes
                         Geometry/Detector)

add_definitions(-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")
add_definitions(-DACTS_GAUDI)

find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_install_headers(JsonWriters)

gaudi_add_module(JsonWriters
                 src/*.cxx
                 INCLUDE_DIRS GaudiKernel ${EIGEN3_HOME} CoreInterfaces Algebra TestInterfaces ParametersBase Surfaces Volumes Detector
                 LINK_LIBRARIES GaudiKernel Surfaces Volumes Detector)
