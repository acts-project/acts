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
        if(NOT ${option})
            string(JOIN " " _acts_option_reason ${ARGN})
            message(
                STATUS
                "Option '${option}' auto-enabled by: ${_acts_option_reason}"
            )
            set_property(
                DIRECTORY
                APPEND
                PROPERTY _acts_auto_enabled_options "${option}"
            )
            set_property(
                DIRECTORY
                PROPERTY _acts_option_reason_${option} "${_acts_option_reason}"
            )
            unset(_acts_option_reason)
        endif()
        # create a regular (directory-scope) variable that should take precedence
        # over a cache entry of the same name. that means that automatically
        # activated options are not stored in the cache.
        set(${option} ON)
    endif()
endmacro()

# Summarize only options changed by set_option_if in this directory. Tracking
# these explicitly keeps internal cache entries out, including on reconfigure.
function(acts_print_auto_enabled_options)
    get_property(_options DIRECTORY PROPERTY _acts_auto_enabled_options)
    if(NOT _options)
        return()
    endif()
    list(REMOVE_DUPLICATES _options)
    list(SORT _options)

    message(STATUS "Automatically enabled options (effective values):")
    foreach(_option IN LISTS _options)
        get_property(_reason DIRECTORY PROPERTY _acts_option_reason_${_option})
        message(STATUS "  ${_option} = ${${_option}} (condition: ${_reason})")
    endforeach()
    message(STATUS "Dependency overrides do not change cached option values.")
endfunction()
