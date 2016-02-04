#ifndef ATS_ATLASPOLICY_H
#define ATS_ATLASPOLICY_H 1

// STL include(s)
#include <cmath>

// ATS include(s)
#include "TrkParameterSet/Parameter_traits.h"
#include "TrkParametersBase/TrackParametersBase.h"
// Trk
#include "TrkEventPrimitives/ParamDefs.h"
#include "TrkSurfaces/Surface.h"

namespace Ats
{
  struct ATLASPolicy
  {
    typedef Trk::ParamDefs par_id_type;
    typedef double par_value_type;
    static constexpr unsigned int N {5};
  };

  template<>
  struct parameter_traits<ATLASPolicy,Trk::loc1>
  {
    typedef local_parameter parameter_type;
  };

  template<>
  struct parameter_traits<ATLASPolicy,Trk::loc2>
  {
    typedef local_parameter parameter_type;
  };

  template<>
  struct parameter_traits<ATLASPolicy,Trk::theta>
  {
    typedef unbound_parameter parameter_type;
  };

  template<>
  struct parameter_traits<ATLASPolicy,Trk::phi>
  {
    static constexpr double pMin(){return -M_PI;}
    static constexpr double pMax(){return M_PI;}
    typedef cyclic_parameter<double,pMin,pMax> parameter_type;
  };

  template<>
  struct parameter_traits<ATLASPolicy,Trk::qOverP>
  {
    typedef unbound_parameter parameter_type;
  };

  template<>
  struct coordinate_transformation<ATLASPolicy>
  {
    typedef AtsVector<typename ATLASPolicy::par_value_type,ATLASPolicy::N> ParVector_t;

    static AtsVectorD<3> parameters2globalPosition(const ParVector_t& pars,const Surface& s)
    {
      AtsVectorD<3> globalPosition;
      s.localToGlobal(AtsVectorD<2>(pars(Trk::loc1),pars(Trk::loc2)),globalPosition,globalPosition);

      return globalPosition;
    }

    static AtsVectorD<3> parameters2globalMomentum(const ParVector_t& pars)
    {
      AtsVectorD<3> momentum;
      double p = fabs(1/pars(Trk::qOverP));
      double phi = pars(Trk::phi);
      double theta = pars(Trk::theta);
      momentum << p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta);

      return momentum;
    }

    static ParVector_t global2curvilinear(const AtsVectorD<3>&,const AtsVectorD<3>& mom,double charge)
    {
      ParVector_t parameters;
      parameters << 0, 0, mom.phi(), mom.theta(), ((fabs(charge) < 1e-4) ? 1. : charge) / mom.mag();

      return parameters;
    }

    static ParVector_t global2parameters(const AtsVectorD<3>& pos,const AtsVectorD<3>& mom,double charge,const Trk::Surface& s)
    {
      AtsVectorD<2> localPosition;
      s.globalToLocal(pos,mom,localPosition);
      ParVector_t result;
      result << localPosition(0),localPosition(1),mom.phi(),mom.theta(),((fabs(charge) < 1e-4) ? 1. : charge) / mom.mag();
      return result;
    }
  };
}  // end of namespace Ats
#endif // ATS_ATLASPOLICY_H
