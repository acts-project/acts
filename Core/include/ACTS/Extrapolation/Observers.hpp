#ifndef ACTS_OBSERVERS_HPP
#define ACTS_OBSERVERS_HPP 1

#include <iostream>
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

struct DebugObserver
{
  template <typename TrackParameters>
  void
  operator()(const TrackParameters& current,
             const TrackParameters& /*previous*/) const
  {
    const auto& pos = current.position();
    std::cout << pos(0) / units::_mm << " " << pos(1) / units::_mm << " "
              << pos(2) / units::_mm << std::endl;
  }
};

}  // namespace Acts

#endif  //  ACTS_OBSERVERS_HPP
