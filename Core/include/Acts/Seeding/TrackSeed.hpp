// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <array>
#include <ostream>
#include <utility>
#include <vector>

#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts {
namespace Seeding {

  /// A set of points and (geometric) parameter estimates.
  ///
  /// @tparam Identifier  A value-like type that is used to link back to the
  ///                     original hits.
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
      , m_phi(Acts::detail::radian_sym(phi))
      , m_theta(theta)
      , m_curvature(curvature)
      , m_ids{{p0.identifier(), ps.identifier()...}}
    {
      static_assert((sizeof...(ps) + 1) == N,
                    "Number of input points must be N");
    }

    const Vector3D&
    position() const
    {
      return m_position;
    }
    double
    phi() const
    {
      return m_phi;
    }
    double
    theta() const
    {
      return m_theta;
    }
    double
    curvature() const
    {
      return m_curvature;
    }
    const std::array<Identifier, N>&
    identifiers() const
    {
      return m_ids;
    }

  private:
    Vector3D m_position;
    double   m_phi, m_theta, m_curvature;
    std::array<Identifier, N> m_ids;
  };

  template <typename Identifier, size_t N>
  std::ostream&
  operator<<(std::ostream& os, const TrackSeed<Identifier, N>& seed)
  {
    using VectorHelpers::phi;
    using VectorHelpers::theta;
    os << "pos=[" << seed.position().x() << ',' << seed.position().y() << ','
       << seed.position().z() << ']';
    os << " phi=" << phi(seed) << " theta=" << theta(seed);
    os << " curvature=" << seed.curvature();
    return os;
  }

  template <typename Identifier, size_t N>
  using TrackSeeds = std::vector<TrackSeed<Identifier, N>>;
  template <typename Identifier>
  using TrackSeeds2 = TrackSeeds<Identifier, 2>;
  template <typename Identifier>
  using TrackSeeds3 = TrackSeeds<Identifier, 3>;
  template <typename Identifier>
  using TrackSeeds4 = TrackSeeds<Identifier, 4>;

}  // namespace Seeding
}  // namespace Acts