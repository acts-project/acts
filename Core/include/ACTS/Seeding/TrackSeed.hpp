// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SEEDING_TRACKSEED_HPP
#define ACTS_SEEDING_TRACKSEED_HPP

#include <array>
#include <ostream>
#include <utility>
#include <vector>

#include "ACTS/Seeding/SpacePoint.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
namespace Seeding {

  /**
   * @brief (Geometric) parameter estimates and a set of points.
   */
  template <typename Identifier, size_t N>
  class TrackSeed
  {
  public:
    static_assert(0 < N, "TrackSeed must have at least one space point.");

    using Point            = SpacePoint<Identifier>;
    using PointIdentifiers = std::array<Identifier, N>;

    template <typename... P>
    TrackSeed(double       phi,
              double       theta,
              double       curvature,
              const Point& p0,
              const P&... ps)
      : m_position(p0.position())
      , m_phi(phi)
      , m_theta(theta)
      , m_curvature(curvature)
      , m_pointIds{p0.identifier(), ps.identifier()...}
    {
      static_assert((sizeof...(ps) + 1) == N,
                    "Number of input points must be N");
    }

    auto
    position() const
    {
      return m_position;
    }
    auto
    phi() const
    {
      return m_phi;
    }
    auto
    theta() const
    {
      return m_theta;
    }
    auto
    curvature() const
    {
      return m_curvature;
    }
    const PointIdentifiers&
    pointIdentifiers() const
    {
      return m_pointIds;
    }

    void
    print(std::ostream& os) const
    {
      os << "x=" << m_position.x() << " y=" << m_position.y()
         << " z=" << m_position.z() << " phi=" << m_phi << " theta=" << m_theta
         << " curvature=" << m_curvature << '\n';
      for (auto id : pointIdentifiers()) os << "  " << id << '\n';
    }

  private:
    Acts::Vector3D   m_position;
    double           m_phi, m_theta, m_curvature;
    PointIdentifiers m_pointIds;
  };

  template <typename Identifier>
  using TrackSeeds2 = std::vector<TrackSeed<Identifier, 2>>;
  template <typename Identifier>
  using TrackSeeds3 = std::vector<TrackSeed<Identifier, 3>>;
  template <typename Identifier>
  using TrackSeeds4 = std::vector<TrackSeed<Identifier, 4>>;

}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_TRACKSEED_HPP
