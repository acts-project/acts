#ifndef ACTS_COORDINATETRANSFORMATIONS_H
#define ACTS_COORDINATETRANSFORMATIONS_H 1

// ACTS includes
#include "Surfaces/Surface.h"
#include "Core/ParameterDefinitions.h"

#ifdef ACTS_COORDINATE_TRANSFORM_PLUGIN

#include ACTS_COORDINATE_TRANSFORM_PLUGIN

#else

namespace Acts
{
  struct coordinate_transformation
  {
    typedef ActsVector<ParValue_t,Acts::NGlobalPars> ParVector_t;

    static ActsVectorD<3> parameters2globalPosition(const ParVector_t& pars,const Surface& s)
    {
      ActsVectorD<3> globalPosition;
      s.localToGlobal(ActsVectorD<2>(pars(Acts::eLOC_1),pars(Acts::eLOC_2)),parameters2globalMomentum(pars),globalPosition);
      return globalPosition;
    }

    static ActsVectorD<3> parameters2globalMomentum(const ParVector_t& pars)
    {
      ActsVectorD<3> momentum;
      double p = fabs(1./pars(Acts::eQOP));
      double phi = pars(Acts::ePHI);
      double theta = pars(Acts::eTHETA);
      momentum << p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta);

      return momentum;
    }

    static ParVector_t global2curvilinear(const ActsVectorD<3>&,const ActsVectorD<3>& mom,double charge)
    {
      ParVector_t parameters;
      parameters << 0, 0, mom.phi(), mom.theta(), ((fabs(charge) < 1e-4) ? 1. : charge) / mom.mag();

      return parameters;
    }

    static ParVector_t global2parameters(const ActsVectorD<3>& pos,const ActsVectorD<3>& mom,double charge,const Surface& s)
    {
      ActsVectorD<2> localPosition;
      s.globalToLocal(pos,mom,localPosition);
      ParVector_t result;
      result << localPosition(0),localPosition(1),mom.phi(),mom.theta(),((fabs(charge) < 1e-4) ? 1. : charge) / mom.mag();
      return result;
    }

    static double parameters2charge(const ParVector_t& pars)
    {
      return (pars(Acts::eQOP) > 0) ? 1. : -1.;
    }
  };
}  // end of namespace Acts
#endif // ACTS_COORDINATE_TRANSFORM_PLUGIN

#endif // ACTS_COORDINATETRANSFORMATIONS_H
