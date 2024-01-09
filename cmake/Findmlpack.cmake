# UNIX paths are standard, no need to specify them.
find_path(mlpack_INCLUDE_DIR
    NAMES mlpack/core.hpp mlpack/prereqs.hpp
    PATHS /opt/mlpack /usr/mlpack
)

find_package_handle_standard_args(mlpack
    REQUIRED_VARS mlpack_INCLUDE_DIR
)

if(mlpack_FOUND)
  set(mlpack_INCLUDE_DIRS ${mlpack_INCLUDE_DIR})
endif()

# Hide internal variables
mark_as_advanced(mlpack_INCLUDE_DIR)
