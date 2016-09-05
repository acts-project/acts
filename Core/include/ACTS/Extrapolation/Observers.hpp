#ifndef ACTS_OBSERVERS_HPP
#define ACTS_OBSERVERS_HPP 1

#include <iostream>
#include <list>
#include <ostream>
#include <tuple>
#include "ACTS/Utilities/Definitions.hpp"
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
    const auto& mom = current.momentum();
    out << pos(0) / units::_mm << " " << pos(1) / units::_mm << " "
        << pos(2) / units::_mm << " " << mom.perp() / units::_GeV << " "
        << mom.phi() << " " << mom.theta() << std::endl;
  }

  mutable std::ostream out{std::cout.rdbuf()};
};

struct PathLengthObserver
{
  typedef struct
  {
    double pathLength = 0.;
  } result_type;

  template <typename TrackParameters>
  void
  operator()(const TrackParameters& current,
             const TrackParameters& previous,
             result_type&           r) const
  {
    const auto& step = current.position() - previous.position();
    r.pathLength += step.norm();
  }
};

struct HitSimulator
{
  struct this_result
  {
    std::list<Vector3D> hits = {};

    void
    add(const this_result& other)
    {
      hits.insert(hits.end(), other.hits.begin(), other.hits.end());
    }

    void
    print(std::ostream& out = std::cout) const
    {
      for (const auto& v : hits)
        out << v(0) / units::_cm << " " << v(1) / units::_cm << " "
            << v(2) / units::_cm << std::endl;
    }
  };

  typedef this_result result_type;

  template <typename TrackParameters>
  void
  operator()(const TrackParameters& current,
             const TrackParameters& previous,
             result_type&           result) const
  {
    double r1 = previous.position().perp();
    double z1 = previous.position().z();
    double r2 = current.position().perp();
    double z2 = current.position().z();

    // check for hit in barrel
    for (const auto& t : barrel) {
      double r    = std::get<0>(t);
      double zmin = std::get<1>(t);
      double zmax = std::get<2>(t);
      if ((r1 - r) * (r2 - r) < 0) {
        double s   = (r - r1) / (r2 - r1);
        auto   pos = s * current.position() + (1 - s) * previous.position();
        if ((zmin - pos.z()) * (zmax - pos.z()) < 0) result.hits.push_back(pos);
      }
    }

    // check for hit in endcaps
    for (const auto& t : endcaps) {
      double z    = std::get<0>(t);
      double rmin = std::get<1>(t);
      double rmax = std::get<2>(t);
      if ((z1 - z) * (z2 - z) < 0) {
        double s   = (z - z1) / (z2 - z1);
        auto   pos = s * current.position() + (1 - s) * previous.position();
        if ((rmin - pos.perp()) * (rmax - pos.perp()) < 0)
          result.hits.push_back(pos);
      }
    }
  }

  std::list<std::tuple<float, float, float>> barrel;
  std::list<std::tuple<float, float, float>> endcaps;
};

}  // namespace Acts

#endif  //  ACTS_OBSERVERS_HPP
