file(MAKE_DIRECTORY "${WORK}")
file(READ "${FIXTURE}" legacy)
# The documentation fixture includes a volume: migration must reject it.
execute_process(
    COMMAND "${MIGRATOR}" "${FIXTURE}" "${WORK}/volume.json"
    RESULT_VARIABLE result
    ERROR_VARIABLE error
)
if(result EQUAL 0 OR NOT error MATCHES "surface material only")
    message(FATAL_ERROR "Expected volume rejection: ${result}: ${error}")
endif()
if(EXISTS "${WORK}/volume.json")
    message(FATAL_ERROR "Failed migration created an output file")
endif()

string(JSON legacy SET "${legacy}" Volumes entries "[]")
string(
    JSON
    keyed_material
    GET "${legacy}"
    Surfaces
    entries
    0
    value
    material
)
string(
    JSON
    legacy
    SET "${legacy}"
    KeyedSurfaces
    "[{\"key\":\"migration-test\",\"geometry_id\":0,\"material\":${keyed_material}}]"
)
file(WRITE "${WORK}/legacy.json" "${legacy}")
execute_process(
    COMMAND "${MIGRATOR}" "${WORK}/legacy.json" "${WORK}/new.json"
    RESULT_VARIABLE result
    ERROR_VARIABLE error
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Migration failed: ${error}")
endif()
file(READ "${WORK}/new.json" migrated)
string(JSON format GET "${migrated}" format)
string(JSON version GET "${migrated}" version)
string(JSON count LENGTH "${migrated}" surfaces)
string(JSON expected_count LENGTH "${legacy}" Surfaces entries)
string(JSON key GET "${migrated}" surfaces ${expected_count} target key)
if(NOT key STREQUAL "migration-test")
    message(FATAL_ERROR "Stable material key was lost")
endif()
math(EXPR expected_count "${expected_count} + 1")
string(
    JSON
    original_thickness
    GET "${legacy}"
    Surfaces
    entries
    0
    value
    material
    data
    0
    0
    thickness
)
string(
    JSON
    migrated_thickness
    GET "${migrated}"
    surfaces
    0
    material
    values
    0
    thickness
)
if(NOT original_thickness STREQUAL migrated_thickness)
    message(FATAL_ERROR "Default migration changed material thickness")
endif()
if(
    NOT format STREQUAL "acts-material-map"
    OR NOT version EQUAL 1
    OR NOT count EQUAL expected_count
)
    message(FATAL_ERROR "Wrong output envelope or surface count")
endif()

execute_process(
    COMMAND
        "${MIGRATOR}" "${WORK}/legacy.json" "${WORK}/quantized.json"
        --material-fraction-bits 16 --compression-level 19 --indentation 2
    RESULT_VARIABLE result
    ERROR_VARIABLE error
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Quantized migration failed: ${error}")
endif()
file(READ "${WORK}/quantized.json" quantized)
string(
    JSON
    quantized_thickness
    GET "${quantized}"
    surfaces
    0
    material
    values
    0
    thickness
)
if(quantized_thickness STREQUAL migrated_thickness)
    message(FATAL_ERROR "Output options had no effect")
endif()

foreach(input IN ITEMS new.json missing.json)
    execute_process(
        COMMAND "${MIGRATOR}" "${WORK}/${input}" "${WORK}/invalid.json"
        RESULT_VARIABLE result
    )
    if(result EQUAL 0)
        message(FATAL_ERROR "Unexpected success for ${input}")
    endif()
endforeach()
foreach(value IN ITEMS 24 -1 abc 16junk)
    execute_process(
        COMMAND
            "${MIGRATOR}" "${WORK}/legacy.json" "${WORK}/invalid.json"
            --material-fraction-bits "${value}"
        RESULT_VARIABLE result
    )
    if(result EQUAL 0)
        message(FATAL_ERROR "Accepted invalid precision ${value}")
    endif()
endforeach()
execute_process(
    COMMAND "${MIGRATOR}" "${WORK}/legacy.json" "${WORK}/legacy.json"
    RESULT_VARIABLE result
)
if(result EQUAL 0)
    message(FATAL_ERROR "Allowed overwriting the input")
endif()
execute_process(COMMAND "${MIGRATOR}" --help RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Help failed")
endif()
