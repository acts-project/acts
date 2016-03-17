################################################################################
# Package: ExtrapolationTest
################################################################################
gaudi_subdir(ExtrapolationTest v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Core/CoreInterfaces
                         Core/Algebra
                         Extrapolation/ExtrapolationInterfaces
                         Extrapolation/ExtrapolationUtils
                         EventData/ParametersBase
                         EventData/TrackParameters
                         EventData/NeutralParameters
                         EventData/EventDataUtils
                         Geometry/Detector
                         Tests/TestInterfaces
                         Geometry/Surfaces)

add_definitions(-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")

find_package(Eigen3 REQUIRED)
find_package(ROOT COMPONENTS MathCore Tree)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_install_headers(ExtrapolationTest)

gaudi_add_module(ExtrapolationTest
                 src/*.cxx
                 INCLUDE_DIRS GaudiKernel CoreInterfaces ${EIGEN3_HOME} ROOT Algebra ExtrapolationInterfaces ExtrapolationUtils ParametersBase TrackParameters NeutralParameters EventDataUtils Detector TestInterfaces Surfaces
                 LINK_LIBRARIES GaudiKernel ROOT ExtrapolationUtils TrackParameters NeutralParameters Detector Surfaces)
