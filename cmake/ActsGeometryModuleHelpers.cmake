function(acts_add_geometry_module target)
    if(NOT DEFINED Acts_GEOMETRY_MODULE_ABI_TAG
       OR Acts_GEOMETRY_MODULE_ABI_TAG STREQUAL ""
    )
        message(
            FATAL_ERROR
            "Acts_GEOMETRY_MODULE_ABI_TAG is not set; cannot configure geometry module target '${target}'."
        )
    endif()

    add_library(${target} SHARED ${ARGN})

    target_link_libraries(${target} PRIVATE Acts::Core)
    target_compile_definitions(
        ${target}
        PRIVATE ACTS_GEOMETRY_MODULE_ABI_TAG="${Acts_GEOMETRY_MODULE_ABI_TAG}"
    )
endfunction()

function(acts_add_dd4hep_geometry_module target)
    if(NOT DEFINED Acts_GEOMETRY_MODULE_ABI_TAG
       OR Acts_GEOMETRY_MODULE_ABI_TAG STREQUAL ""
    )
        message(
            FATAL_ERROR
            "Acts_GEOMETRY_MODULE_ABI_TAG is not set; cannot configure DD4hep geometry module target '${target}'."
        )
    endif()

    add_library(${target} SHARED ${ARGN})

    target_link_libraries(${target} PRIVATE Acts::PluginDD4hep)
    target_compile_definitions(
        ${target}
        PRIVATE ACTS_GEOMETRY_MODULE_ABI_TAG="${Acts_GEOMETRY_MODULE_ABI_TAG}"
    )
endfunction()
