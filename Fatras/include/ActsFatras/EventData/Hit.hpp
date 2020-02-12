// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

namespace ActsFatras {

/// A simulation hit on a surface.
///
/// This is the undigitized, truth hit, i.e. just a recording of the particle
/// state at the surface intersection. Since Fatras is surface-based, the hit
/// position is always constrained to a surface. Depending on the simulated
/// interactions the momentum state before and after might differ and is
/// thus stored as two separate four-vectors.
class Hit {
 public:
  using Scalar = double;
  using Vector4 = Acts::ActsVector<Scalar, 4>;
  using Vector3 = Acts::ActsVector<Scalar, 3>;

  /// Construct default hit with (mostly) invalid information.
  Hit() = default;
  /// Construct from four-position and four-momenta.
  ///
  /// @param geoId      Geometry identifier of the surface
  /// @param particleId Particle identifier of the particle that created the hit
  /// @param pos4       Particle space-time four-vector on the surface
  /// @param before4    Particle four-momentum before the interaction
  /// @param after4     Particle four-momentum after the interaction
  /// @param index_     Hit index along the particle trajectory
  ///
  /// All quantities are given in the global coordinate system. It is the
  /// users responsibility to ensure that the position correspond to a
  /// position on the given surface. The four-vector component order must be
  /// [x,y,z,t] and [px,py,pz,E].
  template <typename Position4, typename Momentum4>
  Hit(Acts::GeometryID geoId, Barcode particleId,
      const Eigen::MatrixBase<Position4>& pos4,
      const Eigen::MatrixBase<Momentum4>& before4,
      const Eigen::MatrixBase<Momentum4>& after4, int32_t index_ = -1)
      : m_geoId(geoId),
        m_particleId(particleId),
        m_index(index_),
        m_pos4(pos4),
        m_before4(before4),
        m_after4(after4) {}
  Hit(const Hit&) = default;
  Hit(Hit&&) = default;
  Hit& operator=(const Hit&) = default;
  Hit& operator=(Hit&&) = default;

  /// Geometry identifier of the hit surface.
  Acts::GeometryID geoId() const { return m_geoId; }
  /// Particle identifier of the particle that generated the hit.
  Barcode particleId() const { return m_particleId; }
  /// Hit index along the particle trajectory.
  ///
  /// @retval negative if the hit index is undefined.
  int32_t index() const { return m_index; }

  /// Space-time position four-vector.
  ///
  /// The component order is [x,y,z,t].
  const Vector4& position4() const { return m_pos4; }
  /// Three-position, i.e. spatial coordinates without the time.
  auto position() const { return m_pos4.head<3>(); }
  /// Time coordinate.
  Scalar time() const { return m_pos4[3]; }

  /// Particle four-momentum before the hit.
  ///
  /// The component order is [px,py,pz,E].
  const Vector4& momentum4Before() const { return m_before4; }
  /// Particle four-momentum after the hit.
  ///
  /// The component order is [px,py,pz,E].
  const Vector4& momentum4After() const { return m_after4; }
  /// Particle direction before the hit.
  Vector3 directionBefore() const { return m_before4.head<3>().normalized(); }
  /// Particle direction after the hit.
  Vector3 directionAfter() const { return m_after4.head<3>().normalized(); }
  /// Average particle direction while passing through the surface.
  Vector3 direction() const {
    auto dir0 = m_before4 / (2 * m_before4.head<3>().norm());
    auto dir1 = m_after4 / (2 * m_after4.head<3>().norm());
    return (dir0 + dir1).head<3>().normalized();
  }
  /// Energy deposited by the hit.
  ///
  /// @retval positive if the particle lost energy when it passed the surface
  /// @retval negative if magic was involved
  Scalar depositedEnergy() const { return m_before4.w() - m_after4.w(); }

 private:
  /// Identifier of the surface.
  Acts::GeometryID m_geoId;
  /// Identifier of the generating particle.
  Barcode m_particleId;
  /// Index of the hit along the particle trajectory.
  int32_t m_index = -1;
  /// Global space-time position four-vector.
  Vector4 m_pos4 = Vector4::Zero();
  /// Global particle energy-momentum four-vector before the hit.
  Vector4 m_before4 = Vector4::Zero();
  /// Global particle energy-momentum four-vector after the hit.
  Vector4 m_after4 = Vector4::Zero();
};

}  // namespace ActsFatras
