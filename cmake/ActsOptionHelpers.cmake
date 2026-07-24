# Provide helper macros to handle option interdependencies.
#
# Ensure an option is set if a condition is met. Can be used to encode
# dependencies between different options, e.g. if OptionA is on, OptionB has
# to be on as well.
#
#     set_option_if(OptionB OptionA)
#
# The macro can take arbitrary conditions as the second parameter, e.g.
#
#     set_option_if(OptionB OptionA AND ConditionB OR ConditionC)
#

macro(set_option_if option)
    if(${ARGN})
        # only announce/record the first time we flip an option from off->on so
        # that user-requested options (already ON in the cache) are not reported
        # as if we had forced them.
        if(NOT ${option})
            string(REPLACE ";" " " _acts_set_option_reason "${ARGN}")
            message(
                STATUS
                "Option '${option}' auto-enabled by: ${_acts_set_option_reason}"
            )
            set_property(
                GLOBAL
                APPEND
                PROPERTY _acts_forced_options "${option}"
            )
            set_property(
                GLOBAL
                PROPERTY
                    _acts_forced_reason_${option} "${_acts_set_option_reason}"
            )
            unset(_acts_set_option_reason)
        endif()
        # create a regular (directory-scope) variable that should take precedence
        # over a cache entry of the same name. that means that automatically
        # activated options are not stored in the cache.
        set(${option} ON)
    endif()
endmacro()

# Snapshot the set of currently-defined ACTS_* cache entries as the canonical
# list of user-facing options. This must be called right after all options have
# been declared, but *before* other machinery (e.g. FetchContent source
# descriptions in ActsExternSources) adds unrelated ACTS_* cache variables, so
# that the summary below only reports genuine options.
function(acts_collect_option_names)
    get_cmake_property(_all_cache_vars CACHE_VARIABLES)
    set(_names "")
    foreach(_var IN LISTS _all_cache_vars)
        if(_var MATCHES "^ACTS_")
            list(APPEND _names "${_var}")
        endif()
    endforeach()
    set_property(GLOBAL PROPERTY _acts_option_names "${_names}")
endfunction()

# Print a consolidated summary of all ACTS_* options and their *effective*
# values at the end of configuration.
#
# This is important because options force-enabled via `set_option_if` are stored
# as directory-scope variables, not cache entries. As a result the cache (and
# tools like `ccmake` / `cmake -LAH`) keep showing the un-forced value, while the
# build actually uses the forced one. This summary reports the value the build
# will really use and flags which options were auto-enabled and why.
function(acts_print_options_summary)
    get_property(_option_names GLOBAL PROPERTY _acts_option_names)
    list(SORT _option_names)
    get_property(_forced GLOBAL PROPERTY _acts_forced_options)

    set(_enabled "")
    set(_disabled "")
    set(_values "")

    foreach(_var IN LISTS _option_names)
        get_property(_type CACHE ${_var} PROPERTY TYPE)
        if(_type STREQUAL "BOOL")
            # `${_var}` reads the regular variable, which shadows the cache entry,
            # i.e. it reflects any auto-enabling done via set_option_if.
            if(${_var})
                if(_forced AND ("${_var}" IN_LIST _forced))
                    get_property(
                        _reason
                        GLOBAL
                        PROPERTY _acts_forced_reason_${_var}
                    )
                    list(
                        APPEND _enabled
                        "  ${_var} = ON   [auto-enabled by: ${_reason}]"
                    )
                else()
                    list(APPEND _enabled "  ${_var} = ON")
                endif()
            else()
                list(APPEND _disabled "${_var}")
            endif()
        elseif(NOT "${${_var}}" STREQUAL "")
            # non-boolean option with a non-default (non-empty) value
            list(APPEND _values "  ${_var} = ${${_var}}")
        endif()
    endforeach()

    list(LENGTH _disabled _n_disabled)
    string(REPLACE ";" ", " _disabled_str "${_disabled}")

    message(STATUS "")
    message(
        STATUS
        "===================== ACTS option summary ====================="
    )
    message(STATUS "Enabled options:")
    foreach(_line IN LISTS _enabled)
        message(STATUS "${_line}")
    endforeach()
    if(_values)
        message(STATUS "Configured values:")
        foreach(_line IN LISTS _values)
            message(STATUS "${_line}")
        endforeach()
    endif()
    message(STATUS "Disabled options (${_n_disabled}): ${_disabled_str}")
    message(
        STATUS
        "  '[auto-enabled by: ...]' = forced on by a dependency; the cache still shows OFF."
    )
    message(
        STATUS
        "==============================================================="
    )
    message(STATUS "")
endfunction()
