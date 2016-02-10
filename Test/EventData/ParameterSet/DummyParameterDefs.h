#ifndef ATS_DUMMYPARAMETERDEFS_H
#define ATS_DUMMYPARAMETERDEFS_H 1

//ATS include(s)
#include "ParameterSet/Parameter_traits.h"

///@cond
namespace Ats
{
  enum class ParDefs {loc1 = 0, loc2 = 1, phi = 2, theta = 3, qop = 4};

  struct ParPolicy
  {
    typedef double    par_value_type;
    typedef ParDefs   par_id_type;
    static constexpr unsigned int N = 5;
  };

  template<>
  struct parameter_traits<ParPolicy,ParDefs::loc1>
  {
    typedef local_parameter parameter_type;
  };

  template<>
  struct parameter_traits<ParPolicy,ParDefs::loc2>
  {
    typedef local_parameter parameter_type;
  };

  template<>
  struct parameter_traits<ParPolicy,ParDefs::theta>
  {
    static constexpr double pMin(){return 0;}
    static constexpr double pMax(){return M_PI;}
    typedef cyclic_parameter<double,pMin,pMax> parameter_type;
  };

  template<>
  struct parameter_traits<ParPolicy,ParDefs::phi>
  {
    static constexpr double pMin(){return -M_PI;}
    static constexpr double pMax(){return M_PI;}
    typedef cyclic_parameter<double,pMin,pMax> parameter_type;
  };

  template<>
  struct parameter_traits<ParPolicy,ParDefs::qop>
  {
    typedef unbound_parameter parameter_type;
  };

}  // end of namespace Ats
/// @endcond

#endif // ATS_DUMMYPARAMETERDEFS_H
