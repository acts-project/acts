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

  /// A set of points and (geometric) parameter estimates.
  ///
  /// @tparam Identifier  A value-like type that is used to link back to the
  ///                     original hits. Must be printable.
  /// @tparam N           Number of hits / space points in this seed.
  template <typename Identifier, size_t N>
  class TrackSeed
  {
  public:
    static_assert(0 < N, "TrackSeed must have at least one space point.");

    template <typename... P>
    TrackSeed(double                        phi,
              double                        theta,
              double                        curvature,
              const SpacePoint<Identifier>& p0,
              const P&... ps)
      : m_position(p0.position())
      , m_phi(phi)
      , m_theta(theta)
      , m_curvature(curvature)
      , m_ids{p0.identifier(), ps.identifier()...}
    {
      static_assert((sizeof...(ps) + 1) == N,
                    "Number of input points must be N");
    }

    const Acts::Vector3D&
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
    const std::array<Identifier, N>&
    identifiers() const
    {
      return m_ids;
    }

    void
    print(std::ostream& os) const
    {
      os << "x=" << m_position.x() << " y=" << m_position.y()
         << " z=" << m_position.z() << " phi=" << m_phi << " theta=" << m_theta
         << " curvature=" << m_curvature << '\n';
      for (const auto& id : m_ids) os << "  " << id << '\n';
    }

  private:
    Acts::Vector3D m_position;
    double         m_phi, m_theta, m_curvature;
    std::array<Identifier, N> m_ids;
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
