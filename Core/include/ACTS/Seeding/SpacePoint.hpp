// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SEEDING_SPACEPOINT_HPP
#define ACTS_SEEDING_SPACEPOINT_HPP

#include <algorithm>
#include <ostream>
#include <vector>

#include "ACTS/Seeding/detail/range.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
namespace Seeding {

  /**
  * @brief 3D space point for track seeding.
  *
  * The identifier is used to link back to the original measurement / hit. It
  * must behave like a value type and can be default initialized.
  */
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

    auto
    position() const
    {
      return m_position;
    }
    auto
    x() const
    {
      return m_position.x();
    }
    auto
    y() const
    {
      return m_position.y();
    }
    auto
    z() const
    {
      return m_position.z();
    }
    auto
    rho() const
    {
      return m_position.head<2>().norm();
    }
    auto
    phi() const
    {
      return m_position.phi();
    }
    auto
    identifier() const
    {
      return m_identifier;
    }

  private:
    Acts::Vector3D m_position;
    Identifier     m_identifier;
  };

  /**
   * @brief A set of space points on a barrel layer with fixed radius.
   *
   * \warning After adding points to the container, `.sort()` must be called to
   * ensure correct ordering. Failure to do so leads to undefined behaviour.
   */
  template <typename Identifier>
  struct BarrelSpacePoints
  {
    std::vector<SpacePoint<Identifier>> points;
    double                              radius;

    /**
     * @brief Subset of points with phi value plus/minus delta around given phi.
     */
    auto
    rangeDeltaPhi(double phi, double delta) const
    {
      return detail::makeRangePhi(points, phi - delta, phi + delta);
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
    os << ' ' << point.identifier();
    return os;
  }

}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_SPACEPOINT_HPP
