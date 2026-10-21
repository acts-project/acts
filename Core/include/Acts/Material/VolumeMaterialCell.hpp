// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/Material.hpp"

#include <functional>
#include <limits>
#include <memory>
#include <utility>

namespace Acts {

/// A region of constant volume material, independent of geometry binning.
///
/// Bounds are lower-inclusive and upper-exclusive in grid coordinates. The
/// optional coordinate transform is shared with the provider so refreshing a
/// cell does not allocate another transform. A homogeneous cell is unbounded.
class VolumeMaterialCell {
 public:
  /// Transform a global position into the coordinates of the material grid.
  using Transform = std::function<Vector3(const Vector3&)>;

  /// Construct an unbounded homogeneous cell.
  /// @param material constant material in the cell
  explicit VolumeMaterialCell(Material material)
      : m_material(std::move(material)) {}

  /// Construct a bounded cell.
  /// @param material constant material in the cell
  /// @param lower lower bounds in grid coordinates
  /// @param upper upper bounds in grid coordinates
  /// @param transform optional global-to-grid coordinate transform
  VolumeMaterialCell(Material material, Vector3 lower, Vector3 upper,
                     std::shared_ptr<const Transform> transform)
      : m_material(std::move(material)),
        m_lower(std::move(lower)),
        m_upper(std::move(upper)),
        m_transform(std::move(transform)) {}

  /// Check whether a global position still lies in this cell.
  /// @param position global lookup position
  /// @return whether the cached material can be reused
  bool isInside(const Vector3& position) const {
    const Vector3 point = m_transform ? (*m_transform)(position) : position;
    return (point.array() >= m_lower.array()).all() &&
           (point.array() < m_upper.array()).all();
  }

  /// @return the constant material in this cell
  const Material& material() const { return m_material; }

 private:
  Material m_material;
  Vector3 m_lower = Vector3::Constant(-std::numeric_limits<double>::infinity());
  Vector3 m_upper = Vector3::Constant(std::numeric_limits<double>::infinity());
  std::shared_ptr<const Transform> m_transform;
};

}  // namespace Acts
