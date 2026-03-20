function(acts_add_geometry_module target)
    if(NOT UNIX)
        message(
            FATAL_ERROR
            "Runtime geometry modules are only supported on Unix-like systems (Linux, macOS). "
            "Cannot configure geometry module target '${target}' on this platform."
        )
    endif()
    if(
        NOT DEFINED ACTS_GEOMETRY_MODULE_ABI_TAG
        OR ACTS_GEOMETRY_MODULE_ABI_TAG STREQUAL ""
    )
        message(
            FATAL_ERROR
            "ACTS_GEOMETRY_MODULE_ABI_TAG is not set; cannot configure geometry module target '${target}'."
        )
    endif()

    add_library(${target} SHARED ${ARGN})

    target_link_libraries(${target} PRIVATE Acts::Core)
    target_compile_definitions(
        ${target}
        PRIVATE ACTS_GEOMETRY_MODULE_ABI_TAG="${ACTS_GEOMETRY_MODULE_ABI_TAG}"
    )
endfunction()
