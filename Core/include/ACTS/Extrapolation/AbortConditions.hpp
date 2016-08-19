#ifndef ACTS_ABORTCONDITIONS_HPP
#define ACTS_ABORTCONDITIONS_HPP 1

#include "ACTS/Extrapolation/Observers.hpp"

namespace Acts {

enum AbortConditions {
  DestinationSurface = 1 << 0,
  MaxMaterial        = 1 << 1,
  MaxPathLength      = 1 << 2
};

template <class Derived, int bitmask>
struct AbortCondition;

template <class Derived>
struct AbortCondition<Derived, DestinationSurface>
{
  typedef struct
  {
  } result_type;

  typedef struct
  {
  } observer_type;

  Derived&
  destinationSurface()
  {
    return static_cast<Derived&>(*this);
  }
};

template <class Derived>
struct AbortCondition<Derived, MaxMaterial>
{
  typedef MaterialObserver           observer_type;
  typedef observer_type::result_type result_type;

  Derived&
  maxMaterial(double m)
  {
    m_maxMaterial = m;
    return static_cast<Derived&>(*this);
  }

  double m_maxMaterial = 0;
};

template <class Derived>
struct AbortCondition<Derived, MaxPathLength>
{
  typedef PathLengthObserver         observer_type;
  typedef observer_type::result_type result_type;

  Derived&
  maxPathLength(double p)
  {
    m_maxPathLength = p;
    return static_cast<Derived&>(*this);
  }

  double m_maxPathLength = 0;
};

}  // namespace Acts
#endif  //  ACTS_ABORTCONDITIONS_HPP
