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
  operator()(const observer_type::result_type& r, TrackParameters&) const
  {
    return (r.pathLength > maxPathLength);
  }
};

}  // namespace Acts
#endif  //  ACTS_ABORTCONDITIONS_HPP
