# Provide helper macros to simplify option interdependencies.
#
# Ensure an option is set if a condition is met. Can be used to encode
# dependencies between different options, e.g. if OptionA is on, OptionB has
# to be on as well.
#
#     set_option_if(OptionB OptionA)
#
# The macro can arbitrary conditions as the second parameter, e.g.
#
#     set_option_if(OptionB OptionA AND ConditionB OR ConditionC)
#

macro(set_option_if option)
  if(${ARGN})
    set(${option} ON)
  endif()
endmacro()
