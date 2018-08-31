// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <algorithm>
#include <ostream>
#include <vector>

#include "Acts/Seeding/detail/cyclic_range.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Seeding {

  /// @brief 3D space point for track seeding.
  ///
  /// @tparam Identifier  A value-like type that is used to link back to the
  ///                     original hits.
  template <typename Identifier>
  class SpacePoint
  {
  public:
    SpacePoint(double x, double y, double z, Identifier id)
      : m_position(x, y, z), m_identifier(id)
    {
    }
    SpacePoint(const Acts::Vector3D& xyz, Identifier id)
      : m_position(xyz), m_identifier(id)
    {
    }

    const Vector3D&
    position() const
    {
      return m_position;
    }
    double
    x() const
    {
      return m_position.x();
    }
    double
    y() const
    {
      return m_position.y();
    }
    double
    z() const
    {
      return m_position.z();
    }
    double
    rho() const
    {
      return m_position.head<2>().norm();
    }
    double
    phi() const
    {
      return VectorHelpers::phi(m_position);
    }
    const Identifier&
    identifier() const
    {
      return m_identifier;
    }

  private:
    Acts::Vector3D m_position;
    Identifier     m_identifier;
  };

  /// @brief A set of space points on a barrel layer with fixed radius.
  ///
  /// @warning After adding points to the container, `.sort()` must be called
  ///          to ensure correct ordering. Failure to do so leads to undefined
  ///          behaviour.
  template <typename Identifier>
  struct BarrelSpacePoints
  {
    std::vector<SpacePoint<Identifier>> points;
    double                              radius;

    /// Subset of points with phi value plus/minus delta around given phi.
    auto
    rangeDeltaPhi(double phi, double delta) const
    {
      return Acts::detail::makeRangePhi(points, phi - delta, phi + delta);
    }

    void
    sort()
    {
      auto compare
          = [](const auto p0, const auto& p1) { return (p0.phi() < p1.phi()); };
      std::sort(points.begin(), points.end(), compare);
    }

    BarrelSpacePoints(double radius_) : radius(radius_) {}
  };

  template <typename Identifier>
  inline std::ostream&
  operator<<(std::ostream& os, const SpacePoint<Identifier>& point)
  {
    os << "rho=" << point.rho() << " phi=" << point.phi() << " z=" << point.z();
    return os;
  }

}  // namespace Seeding
}  // namespace Acts