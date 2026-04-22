# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Let the user know what's happening.
message(STATUS "Building 'detray::io' component")

# Set up the I/O utils library.
file(
    GLOB _detray_io_utils_public_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/io/utils/*.hpp"
)
detray_add_library( detray_io_utils io_utils
   ${_detray_io_utils_public_headers}
)

# Set up the core I/O library.
file(
    GLOB _detray_io_public_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/io/backend/*.hpp"
    "include/detray/io/backend/detail/*.hpp"
    "include/detray/io/frontend/*.hpp"
    "include/detray/io/frontend/detail/*.hpp"
    "include/detray/io/frontend/impl/*.hpp"
    "include/detray/io/json/*.hpp"
    "include/detray/io/json/detail/*.hpp"
)
detray_add_library( detray_io io
   ${_detray_io_public_headers}
)
target_link_libraries(
    detray_io
    INTERFACE
        nlohmann_json::nlohmann_json
        vecmem::core
        detray::core
        detray::io_utils
)

# Set up libraries using particular algebra plugins.
detray_add_library( detray_io_array io_array )
target_link_libraries(
    detray_io_array
    INTERFACE detray::io detray::algebra_array
)

if(DETRAY_EIGEN_PLUGIN)
    detray_add_library( detray_io_eigen io_eigen )
    target_link_libraries(
        detray_io_eigen
        INTERFACE detray::io detray::algebra_eigen
    )
endif()

if(DETRAY_SMATRIX_PLUGIN)
    detray_add_library( detray_io_smatrix io_smatrix )
    target_link_libraries(
        detray_io_smatrix
        INTERFACE detray::io detray::algebra_smatrix
    )
endif()

if(DETRAY_VC_AOS_PLUGIN)
    detray_add_library( detray_io_vc_aos io_vc_aos )
    target_link_libraries(
        detray_io_vc_aos
        INTERFACE detray::io detray::algebra_vc_aos
    )
endif()

# Set up the csv data I/O library.
file(
    GLOB _detray_csv_io_public_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/io/csv/*.hpp"
)
detray_add_library( detray_csv_io csv_io
   ${_detray_csv_io_public_headers}
)
target_link_libraries(detray_csv_io INTERFACE detray::io_utils)

# Test the public headers of the detray I/O libraries.
if(BUILD_TESTING AND DETRAY_BUILD_TESTING)
    string(
        REPLACE "include/"
        ""
        _detray_io_utils_public_headers
        "${_detray_io_utils_public_headers}"
    )

    string(
        REPLACE "include/"
        ""
        _detray_io_public_headers
        "${_detray_io_public_headers}"
    )
    detray_test_public_headers(detray_io_array
      ${_detray_io_public_headers}
    )
    if(DETRAY_EIGEN_PLUGIN)
        detray_test_public_headers(detray_io_eigen
         ${_detray_io_public_headers}
        )
    endif()
    if(DETRAY_SMATRIX_PLUGIN)
        detray_test_public_headers(detray_io_smatrix
         ${_detray_io_public_headers}
        )
    endif()
    if(DETRAY_VC_AOS_PLUGIN)
        detray_test_public_headers(detray_io_vc_aos
         ${_detray_io_public_headers}
        )
    endif()

    string(
        REPLACE "include/"
        ""
        _detray_csv_io_public_headers
        "${_detray_csv_io_public_headers}"
    )
endif()

# Clean up.
unset(_detray_io_public_headers)
