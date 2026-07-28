// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/TransformRange.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

namespace ActsAlignment {

/// @brief A logical group of surfaces that share alignment degrees of freedom.
///
/// AlignableStructure aggregates a set of modules (surfaces) into a single
/// alignable unit parameterised by the 6 standard rigid-body DoFs
/// (Center0..2, Rotation0..2) plus optional per-DoF rigidity constraints.
/// It owns no transform of its own; alignment deltas live in the
/// AlignmentContext store, exactly like the per-module case.
class AlignableStructure {
 public:
  /// @brief Constructor
  /// @param id A unique identifier for this structure. Must be distinct from
  ///           every other structure in the same hierarchy and must not collide
  ///           with any surface @c GeometryIdentifier already present in the
  ///           alignment store (surface IDs always carry non-zero sensitive
  ///           bits, so keeping the sensitive field zero here is sufficient).
  ///           The @c GeometryIdentifier type is chosen because it is the key
  ///           type of @c GeoIdAlignmentStore; the value need not correspond to
  ///           any existing geometry object.
  /// @todo Define a project-wide strategy for assigning structure IDs.
  explicit AlignableStructure(Acts::GeometryIdentifier id) : m_id(id) {}

  /// Disallow copy construction to ensure identity uniqueness
  AlignableStructure(const AlignableStructure&) = delete;
  AlignableStructure& operator=(const AlignableStructure&) = delete;

  /// Type alias for a const dereferencing range over associated surfaces
  using SurfaceRange = Acts::detail::TransformRange<
      Acts::detail::ConstDereference,
      const std::vector<std::shared_ptr<Acts::Surface>>>;

  /// @brief Add a surface to this alignable structure
  /// @param surface The surface to add
  void addSurface(std::shared_ptr<Acts::Surface> surface) {
    m_surfaces.push_back(std::move(surface));
  }

  /// @brief Add a child structure for nested hierarchies
  /// @param child Shared pointer to a child structure
  void addChild(std::shared_ptr<AlignableStructure> child) {
    m_children.push_back(std::move(child));
  }

  /// @brief Get the geometric identifier
  /// @return The ID of the structure
  Acts::GeometryIdentifier geometryId() const { return m_id; }

  /// @brief Access the surfaces directly associated with this structure.
  /// @return A dereferencing range yielding @c const Acts::Surface&
  SurfaceRange surfaces() const { return SurfaceRange{m_surfaces}; }

  /// @brief Access the list of child structures
  /// @return A vector of shared pointers to child structures
  const std::vector<std::shared_ptr<AlignableStructure>>& children() const {
    return m_children;
  }

  /// @brief Access the constraints map (rigidity)
  ///
  /// The map key is the standard alignment index. The value is the variance
  /// constraint \f$\sigma_i^2\f$ contributing \f$W_{ii} = 1/\sigma_i^2\f$ to
  /// the Hessian diagonal. Setting \f$\sigma_i \to 0\f$ locks DoF @c i.
  /// A DoF that is floating (set in @c alignmentMask) but absent from this
  /// map is unconstrained — no regularization term is added for it.
  /// @return Reference to the constraints map
  std::unordered_map<Acts::AlignmentIndices, double>& constraints() {
    return m_constraints;
  }
  /// @brief Access the constraints map (const overload)
  /// @return Const reference to the constraints map
  const std::unordered_map<Acts::AlignmentIndices, double>& constraints()
      const {
    return m_constraints;
  }

  /// @brief Access the alignment mask (floating DoFs)
  /// @return Reference to the mask
  AlignmentMask& alignmentMask() { return m_mask; }
  /// @brief Access the alignment mask (const overload)
  /// @return The mask value
  AlignmentMask alignmentMask() const { return m_mask; }

 private:
  /// The unique identifier for this structure
  Acts::GeometryIdentifier m_id;

  /// The collection of sensors governed by this structure
  std::vector<std::shared_ptr<Acts::Surface>> m_surfaces;

  /// The child structures for nested hierarchies
  std::vector<std::shared_ptr<AlignableStructure>> m_children;

  /// Map of parameter index -> variance constraint (rigidity)
  std::unordered_map<Acts::AlignmentIndices, double> m_constraints;

  /// Bitmask for determining which of the 6 DoFs are floating
  AlignmentMask m_mask{AlignmentMask::None};
};

}  // namespace ActsAlignment
