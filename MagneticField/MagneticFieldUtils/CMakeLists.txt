################################################################################
# Package: MagneticFieldUtils
################################################################################
gaudi_subdir(MagneticFieldUtils v1r0)


gaudi_depends_on_subdirs(Core/Algebra)

find_package(Boost)
find_package(Eigen3 REQUIRED)

SET( EIGEN3_HOME ${EIGEN_ROOT_DIR}/${EIGEN_INCLUDE_DIR} )

gaudi_add_library(MagneticFieldUtils
		  src/*.cxx
		  PUBLIC_HEADERS MagneticFieldUtils
		  INCLUDE_DIRS ${EIGEN3_HOME} Algebra Boost
		  LINK_LIBRARIES Boost)
