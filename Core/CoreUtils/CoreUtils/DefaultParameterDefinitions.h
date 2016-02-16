#ifndef ATS_DEFAULTPARAMETERDEFINITIONS_H
#define ATS_DEFAULTPARAMETERDEFINITIONS_H 1

// STL include(s)
#include <cmath>

// ATS include(s)
#include "CoreUtils/ParameterTypes.h"

namespace Ats
{
  enum ParDef : unsigned int
  {
    eLOC_1   = 0, ///< first coordinate in local surface frame
    eLOC_2   = 1, ///< second coordinate in local surface frame
    eLOC_R = eLOC_1,
    eLOC_PHI = eLOC_2,
    eLOC_RPHI = eLOC_1,
    eLOC_Z = eLOC_2,
    eLOC_X = eLOC_1,
    eLOC_Y = eLOC_2,
    eLOC_D0 = eLOC_1,
    eLOC_Z0 = eLOC_2,
    ePHI    = 2, ///< phi direction of momentum in global frame
    eTHETA  = 3, ///< theta direction of momentum in global frame
    eQOP    = 4,  ///< charge/momentum for charged tracks, for neutral tracks it is 1/momentum
    NGlobalPars
  };

  typedef ParDef ParID_t;
  typedef double ParValue_t;

  template<ParID_t>
  struct par_type;

  template<ParID_t par>
  using par_type_t = typename par_type<par>::type;

  template<>
  struct par_type<ParDef::eLOC_1>
  {
    typedef local_parameter type;
  };

  template<>
  struct par_type<ParDef::eLOC_2>
  {
    typedef local_parameter type;
  };

  template<>
  struct par_type<ParDef::ePHI>
  {
    static constexpr double pMin(){return -M_PI;}
    static constexpr double pMax(){return M_PI;}
    typedef cyclic_parameter<double,pMin,pMax> type;
  };

  template<>
  struct par_type<ParDef::eTHETA>
  {
    static constexpr double pMin(){return 0;}
    static constexpr double pMax(){return M_PI;}
    typedef bound_parameter<double,pMin,pMax> type;
  };

  template<>
  struct par_type<ParDef::eQOP>
  {
    typedef unbound_parameter type;
  };
}  // end of namespace Ats

#endif // ATS_DEFAULTPARAMETERDEFINITIONS_H
