function(acts_add_dd4hep_geometry_module target)
    if(NOT UNIX)
        message(
            FATAL_ERROR
            "Runtime geometry modules are only supported on Unix-like systems (Linux, macOS). "
            "Cannot configure DD4hep geometry module target '${target}' on this platform."
        )
    endif()
    if(
        NOT DEFINED ACTS_GEOMETRY_MODULE_ABI_TAG
        OR ACTS_GEOMETRY_MODULE_ABI_TAG STREQUAL ""
    )
        message(
            FATAL_ERROR
            "ACTS_GEOMETRY_MODULE_ABI_TAG is not set; cannot configure DD4hep geometry module target '${target}'."
        )
    endif()
    if(NOT DEFINED DD4hep_VERSION OR DD4hep_VERSION STREQUAL "")
        message(
            FATAL_ERROR
            "DD4hep_VERSION is not set; cannot configure DD4hep geometry module target '${target}'."
        )
    endif()

    set(_acts_dd4hep_geometry_module_abi_tag
        "${ACTS_GEOMETRY_MODULE_ABI_TAG}|dd4hep-${DD4hep_VERSION}"
    )

    add_library(${target} SHARED ${ARGN})

    target_link_libraries(${target} PRIVATE Acts::PluginDD4hep)
    target_compile_definitions(
        ${target}
        PRIVATE
            ACTS_GEOMETRY_MODULE_ABI_TAG="${_acts_dd4hep_geometry_module_abi_tag}"
    )
endfunction()
