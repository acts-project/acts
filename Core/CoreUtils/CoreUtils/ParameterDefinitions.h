#ifndef ATS_PARAMETERDEFINITIONS_H
#define ATS_PARAMETERDEFINITIONS_H 1

#ifdef ATS_PARAMETER_DEFINITIONS_PLUGIN
#include ATS_PARAMETER_DEFINITIONS_PLUGIN
#endif

// testing requirements on parameter definitions
#include <type_traits>

// typedefs for parameter identifier and parameter value must be present
static_assert(std::is_enum<Ats::ParID_t>::value,"'ParID_t' is not an enum type");
static_assert(std::is_floating_point<Ats::ParValue_t>::value,"'ParValue_t' is not floating point type");

// parameter ID type must be convertible to unsigned int
static_assert(std::is_convertible<Ats::ParID_t,unsigned int>::value,"'ParID_t' is not convertible to unsigned int");

// number of global parameter must be at least 2 (for the two local parameters)
static_assert(Ats::NGlobalPars > 1,"total number of global parameters must be >= 2");

// several constants for the local parameters need to be defined
static_assert(Ats::eLOC_1 != Ats::eLOC_2,"local parameters must have different IDs");
static_assert(Ats::eLOC_R == Ats::eLOC_1 or Ats::eLOC_R == Ats::eLOC_2,"local radius must be a local parameter");
static_assert(Ats::eLOC_PHI == Ats::eLOC_1 or Ats::eLOC_PHI == Ats::eLOC_2,"local phi must be a local parameter");
static_assert(Ats::eLOC_RPHI == Ats::eLOC_1 or Ats::eLOC_RPHI == Ats::eLOC_2,"local r x phi must be a local parameter");
static_assert(Ats::eLOC_Z == Ats::eLOC_1 or Ats::eLOC_Z == Ats::eLOC_2,"local z must be a local parameter");
static_assert(Ats::eLOC_X == Ats::eLOC_1 or Ats::eLOC_X == Ats::eLOC_2,"local x must be a local parameter");
static_assert(Ats::eLOC_Y == Ats::eLOC_1 or Ats::eLOC_Y == Ats::eLOC_2,"local y must be a local parameter");
static_assert(Ats::eLOC_D0 == Ats::eLOC_1 or Ats::eLOC_D0 == Ats::eLOC_2,"d0 must be a local parameter");
static_assert(Ats::eLOC_Z0 == Ats::eLOC_1 or Ats::eLOC_Z0 == Ats::eLOC_2,"z0 must be a local parameter");

// check for par_type_t definition
static_assert(sizeof(Ats::par_type_t<Ats::eLOC_1>) > 0,"'par_type_t' is not defined");

#endif //  ATS_PARAMETERDEFINITIONS_H
