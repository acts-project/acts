################################################################################
# Package: FatrasTools
################################################################################
gaudi_subdir(FatrasTools v1r0)

gaudi_depends_on_subdirs(GaudiKernel
                         Core/CoreInterfaces
                         Core/Algebra
                         Extensions/Fatras/FatrasInterfaces
                         EventData/EventDataUtils
                         Extrapolation/ExtrapolationUtils
                         Geometry/Material
                         Extensions/Fatras/FatrasUtils)

add_definitions(-DACTS_PARAMETER_DEFINITIONS_PLUGIN=\"CoreUtils/DefaultParameterDefinitions.h\")
add_definitions(-DACTS_GAUDI)

find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_install_headers(FatrasTools)

gaudi_add_module(FatrasTools
                src/*.cxx
                INCLUDE_DIRS GaudiKernel CoreInterfaces Algebra ${EIGEN3_HOME} FatrasInterfaces ExtrapolationUtils Material FatrasUtils
                LINK_LIBRARIES GaudiKernel ExtrapolationUtils Material FatrasUtils)


