#ifndef ACTS_ABORTCONDITIONS_HPP
#define ACTS_ABORTCONDITIONS_HPP 1

#include "ACTS/Extrapolation/Observers.hpp"

namespace Acts {

class Surface;

struct DestinationSurface
{
  void
  setTargetSurface(const Surface& target)
  {
    m_pTarget = &target;
  }

private:
  const Surface* m_pTarget = 0;
};

// struct MaxMaterial
//{
//  typedef MaterialObserver           observer_type;
//  typedef observer_type::result_type result_type;
//
//  double maxMaterial = 0;
//};
//
struct MaxPathLength
{
  typedef PathLengthObserver observer_type;

  double maxPathLength = 0;

  template <typename TrackParameters>
  bool
  operator()(const observer_type::result_type& r,
             TrackParameters&,
             double& stepMax) const
  {
    // adjust maximum step size
    stepMax = maxPathLength - r.pathLength;

    return (r.pathLength >= maxPathLength);
  }
};

struct MaxRadius
{
  double maxRadius = 0;

  template <typename TrackParameters>
  bool
  operator()(const TrackParameters& pars, double& stepMax) const
  {
    // adjust maximum step size
    // dR/dS = pT/p -> dS = dR * p / pT = dR * sqrt( 1 + pZ^2 / pT^2)
    const double deltaR       = pars.position().perp() - maxRadius;
    const double pT2          = pars.momentum().perp2();
    const double pZ2          = pars.momentum().z() * pars.momentum().z();
    const double inv_sintheta = sqrt(1 + pZ2 / pT2);

    stepMax = -deltaR * inv_sintheta;

    return (deltaR > 0);
  }
};

}  // namespace Acts
#endif  //  ACTS_ABORTCONDITIONS_HPP
