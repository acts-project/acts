// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/TransformRange.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <functional>
#include <memory>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ActsAlignment {

/// @brief A logical group of surfaces that share alignment degrees of freedom.
///
/// AlignableStructure aggregates a set of modules (surfaces) into a single
/// alignable unit parameterised by the 6 standard rigid-body DoFs
/// (Center0..2, Rotation0..2) plus optional per-DoF rigidity constraints.
/// The structure frame defines its origin (rotation pivot) and local axes
/// in global coordinates. It is supplied either as a fixed transform or
/// as a TransformProvider callable evaluated with a GeometryContext.
/// The returned frame must be consistent with the member-surface transforms
/// evaluated using the same GeometryContext.
/// Both structure translations and rotations use these local axes, following
/// @c Acts::detail::compositeToComponentJacobian.
/// Querying the frame does not move any member surface.
/// A structure constructed with only an ID is a frame-less declaration of
/// membership, suitable for hierarchy validation. Numerical alignment also
/// requires a frame, which can be supplied later with setTransform().
class AlignableStructure {
 public:
  /// Resolve the current structure-to-global rigid transform for a context.
  /// The callable must be safe for concurrent reads and must not mutate shared
  /// geometry. Captured conditions data must outlive its use. Returning by
  /// value avoids exposing references into conditions caches.
  using TransformProvider =
      std::function<Acts::Transform3(const Acts::GeometryContext&)>;

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

  /// Construct a structure with a context-independent frame.
  /// This fixes the frame returned by localToGlobalTransform(), not the
  /// alignment degrees of freedom, which are controlled by alignmentMask().
  /// @param id Identifier with the same uniqueness requirements as the ID-only constructor
  /// @param transform Rigid structure-to-global transform, copied into the provider
  /// @note This frame is independent of the context. Use the constructor taking
  /// a callable when the structure's pose depends on an IOV or alignment state.
  /// The caller chooses the pivot; it is not inferred from member surface
  /// centers. A solver using this as an initial frame must maintain subsequent
  /// frame updates separately; the supplied frame is used until replaced with
  /// setTransform().
  AlignableStructure(Acts::GeometryIdentifier id,
                     const Acts::Transform3& transform)
      : AlignableStructure(id, [transform](const Acts::GeometryContext&) {
          return transform;
        }) {}

  /// Construct a structure with a context-dependent frame.
  /// @param id Identifier with the same uniqueness requirements as the ID-only constructor
  /// @param transformProvider Callable returning the structure-to-global rigid transform
  /// @throws std::invalid_argument if the provider is empty
  /// @note The returned child frame must map child-local to global coordinates.
  /// The callable may compose a parent-to-global frame with a child-to-parent
  /// transform, both evaluated for the same context. Hierarchy membership does
  /// not perform this composition or update frames automatically.
  AlignableStructure(Acts::GeometryIdentifier id,
                     TransformProvider transformProvider)
      : m_id(id) {
    setTransform(std::move(transformProvider));
  }

  /// Add or replace the frame with a context-independent transform.
  /// @param transform Rigid structure-to-global transform, copied into the provider
  /// @note This configures the frame lookup; it does not move member surfaces.
  /// The caller must ensure there are no concurrent readers or writers.
  void setTransform(const Acts::Transform3& transform) {
    setTransform(
        [transform](const Acts::GeometryContext&) { return transform; });
  }

  /// Add or replace the context-dependent frame provider.
  /// @param transformProvider Callable returning the structure-to-global rigid transform
  /// @throws std::invalid_argument if the provider is empty; the existing provider is preserved
  /// @note This configures the frame lookup; it does not move member surfaces.
  /// The caller must ensure there are no concurrent readers or writers. Use
  /// GeometryContext with the installed provider for ordinary IOV changes.
  void setTransform(TransformProvider transformProvider) {
    if (!transformProvider) {
      throw std::invalid_argument("Structure transform provider is empty");
    }
    m_transformProvider = std::move(transformProvider);
  }

  /// Whether a frame provider is currently installed.
  /// @return True for a fixed or contextual frame, false if no frame was supplied
  bool hasTransform() const { return static_cast<bool>(m_transformProvider); }

  /// Resolve this structure's origin and axes in the global coordinate system.
  /// @param gctx The same geometry context used to evaluate member surfaces
  /// @return A copy of the current rigid structure-to-global transform
  /// @throws std::logic_error if this structure has no frame provider
  /// @note Provider exceptions propagate to the caller. The provider is queried
  /// on every call; querying does not cache or modify any transform.
  /// Corrections must be applied consistently to the structure frame and its
  /// member surfaces by the owner of the alignment state before evaluating the
  /// next iteration.
  Acts::Transform3 localToGlobalTransform(
      const Acts::GeometryContext& gctx) const {
    if (!m_transformProvider) {
      throw std::logic_error("Alignable structure has no frame provider");
    }
    return m_transformProvider(gctx);
  }

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
  /// @return A read-only view of the shared pointers to child structures
  std::span<const std::shared_ptr<AlignableStructure>> children() const {
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

  /// Frame lookup configured by the caller; contextual state belongs to the
  /// caller
  TransformProvider m_transformProvider;

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
