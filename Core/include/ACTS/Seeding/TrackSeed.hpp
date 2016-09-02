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

namespace Acts {
namespace Seeding {

  /**
   * @brief (Geometric) parameter estimates and a set of space points.
   *
   * The parameter estimates are always defined at the first space point.
   */
  template <typename Identifier, size_t N>
  class TrackSeed
  {
  public:
    static_assert(0 < N, "TrackSeed must have at least one space point.");

    static constexpr size_t kNumPoints = N;

    template <typename... P>
    TrackSeed(double phi, double theta, double curvature, P&&... ps)
      : m_phi(phi)
      , m_theta(theta)
      , m_curvature(curvature)
      , m_points{std::forward<P>(ps)...}
    {
      static_assert(sizeof...(ps) == N, "Number of input points must be N");
    }
    template <typename... P>
    TrackSeed(const Vector3D direction, double curvature, P&&... ps)
      : TrackSeed(direction.phi(),
                  direction.theta(),
                  curvature,
                  std::forward<P>(ps)...)
    {
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
    template <size_t INDEX>
    constexpr const SpacePoint<Identifier>&
    point() const
    {
      return std::get<INDEX>(m_points);
    }
    const SpacePoint<Identifier>&
    point(size_t index) const
    {
      return m_points[index];
    }

    void print(std::ostream& os) const {
      os << "phi=" << m_phi << " theta=" << m_theta
         << " curvature=" << m_curvature << '\n';
      for (const auto& point : m_points)
        os << "  " << point << '\n';
    }

  private:
    double m_phi, m_theta, m_curvature;
    std::array<const SpacePoint<Identifier>, N> m_points;
  };

  template <typename Identifier>
  using TrackSeeds2 = std::vector<TrackSeed<Identifier, 2>>;
  template <typename Identifier>
  using TrackSeeds3 = std::vector<TrackSeed<Identifier, 3>>;
  template <typename Identifier>
  using TrackSeeds4 = std::vector<TrackSeed<Identifier, 4>>;

  template <typename Identifier, size_t N>
  inline std::ostream&
  operator<<(std::ostream& os, const TrackSeed<Identifier, N>& seed)
  {
    os << "phi=" << seed.phi() << " theta=" << seed.theta()
       << " kappa=" << seed.curvature() << " n_points=" << N;
    return os;
  }

}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_TRACKSEED_HPP
