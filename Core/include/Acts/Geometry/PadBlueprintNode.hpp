// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <optional>

namespace Acts {

/// Wraps a single child blueprint node and pads it into a larger volume whose
/// dimensions are evaluated at construction time using the child's extent plus
/// a configurable envelope.
/// @note This node can only have a single child. This is not an error during
///       tree building, but during geometry construction.
/// It defers most of the functionality to @ref Acts::StaticBlueprintNode,
/// and only implements the build phase to perform the padding.
///
/// Two orthogonal, opt-in settings control where the padded volume ends up:
/// - @ref setReferenceAxis picks the @b frame the enclosure is built in. By
///   default it is the child's own frame. When a reference axis is set, a
///   cylinder child's transverse offset from that axis is absorbed into radial
///   growth and the enclosure is re-centered onto the axis, so a displaced
///   cylindrical subtree can be dropped into a strictly co-axial
///   @ref Acts::CylinderVolumeStack without weakening its alignment checks.
///   Axis alignment is cylinder-only.
/// - @ref setCentering picks how the volume is centered along the directions
///   its bounds express symmetrically (z for a cylinder, x/y/z for a cuboid).
///   The default @ref Centering::Centered keeps the volume centered on the
///   child and requires a symmetric envelope there; @ref Centering::FitBounds
///   allows an asymmetric envelope and shifts the enclosure so the child stays
///   contained.
class PadBlueprintNode final : public StaticBlueprintNode {
 public:
  /// Controls how the padded volume is centered along the directions its bounds
  /// express symmetrically (z for a cylinder, x/y/z for a cuboid).
  enum class Centering {
    /// Keep the padded volume centered on the child (default). The envelope
    /// must be symmetric in those directions, otherwise @ref build throws: the
    /// volume cannot move to absorb the asymmetry.
    Centered,
    /// Place the padded volume at the midpoint of the *expanded* extent, so an
    /// asymmetric envelope simply shifts the enclosure off the child while
    /// keeping it fully contained.
    FitBounds,
  };

  /// Main constructor for the padding node.
  /// @param name The name of the padded volume.
  /// @param envelope The envelope to apply to the child node's extent to create the padded volume.
  explicit PadBlueprintNode(
      const std::string& name,
      const ExtentEnvelope& envelope = ExtentEnvelope::Zero());

  ~PadBlueprintNode() override = default;

  /// @copydoc BlueprintNode::build
  /// Build-phase of the blueprint construction. Returns the padded volume.
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// Set the padding envelope applied to the child's extent.
  /// @param envelope The per-side envelope to apply
  /// @return Reference to this node for chaining
  PadBlueprintNode& setEnvelope(const ExtentEnvelope& envelope);

  /// Get the padding envelope.
  /// @return The configured envelope
  const ExtentEnvelope& envelope() const;

  /// Align the enclosure to an external reference axis instead of the child's
  /// own frame (cylinder children only). The child's transverse offset from
  /// @p axis is absorbed into radial growth and the enclosure is re-centered
  /// onto @p axis; the child keeps its displaced placement. This transverse
  /// re-centering happens irrespective of the @ref Centering setting, which
  /// still governs the axial (z) direction.
  /// @param axis Transform of the reference axis frame
  /// @return Reference to this node for chaining
  PadBlueprintNode& setReferenceAxis(const Transform3& axis);

  /// Drop the reference axis and pad in the child's own frame (the default).
  /// @return Reference to this node for chaining
  PadBlueprintNode& clearReferenceAxis();

  /// Get the reference axis, if one is configured.
  /// @return The optional reference axis transform
  const std::optional<Transform3>& referenceAxis() const;

  /// Set how the volume is centered along its symmetric directions.
  /// @param centering The centering mode
  /// @return Reference to this node for chaining
  PadBlueprintNode& setCentering(Centering centering);

  /// Get the centering mode.
  /// @return The configured centering mode
  Centering centering() const;

  /// Create a volume that encloses @p inner, enlarged by @p envelope.
  ///
  /// Without a @p referenceAxis the enclosure is built in the child's own
  /// frame; with one (cylinder children only) it is re-centered onto that axis
  /// and its radial bounds are grown by the child's radial offset. The
  /// @p centering mode governs the symmetric directions (z for cylinders,
  /// x/y/z for cuboids): @ref Centering::Centered keeps the volume centered on
  /// the child and requires a symmetric envelope there, while
  /// @ref Centering::FitBounds allows an asymmetric envelope and shifts the
  /// enclosure to the midpoint of the expanded extent.
  /// @param gctx The geometry context
  /// @param inner The volume to enclose
  /// @param envelope The envelope to add to the bounds of @p inner
  /// @param name The name of the padded volume
  /// @param referenceAxis Optional reference axis for axis-aligned enclosure
  /// @param centering Centering mode for the symmetric directions
  /// @param logger The logger to use
  /// @return The padded volume enclosing @p inner
  /// @throws std::logic_error if @p inner has unsupported bounds, if a
  ///         @ref Centering::Centered envelope is asymmetric where it must not
  ///         be, or if a reference axis is used with a non-cylinder or tilted
  ///         child
  static std::unique_ptr<TrackingVolume> padded(
      const GeometryContext& gctx, const Volume& inner,
      const ExtentEnvelope& envelope, const std::string& name,
      const std::optional<Transform3>& referenceAxis = std::nullopt,
      Centering centering = Centering::Centered,
      const Logger& logger = Acts::getDummyLogger());

  /// Overload keeping the pre-reference-axis call signature source compatible,
  /// i.e. padding in the child's own frame with @ref Centering::Centered.
  /// @param gctx The geometry context
  /// @param inner The volume to enclose
  /// @param envelope The envelope to add to the bounds of @p inner
  /// @param name The name of the padded volume
  /// @param logger The logger to use
  /// @return The padded volume enclosing @p inner
  /// @throws std::logic_error if @p inner has unsupported bounds, or if the
  ///         envelope is asymmetric where it must not be
  static std::unique_ptr<TrackingVolume> padded(const GeometryContext& gctx,
                                                const Volume& inner,
                                                const ExtentEnvelope& envelope,
                                                const std::string& name,
                                                const Logger& logger);

 private:
  ExtentEnvelope m_envelope;
  std::string m_name;
  std::optional<Transform3> m_referenceAxis;
  Centering m_centering = Centering::Centered;
};

namespace Experimental {
/// @deprecated The blueprint geometry moved out of the `Acts::Experimental`
///             namespace. Use @ref Acts::PadBlueprintNode instead. This alias
///             is kept for backward compatibility and will be removed.
using PadBlueprintNode [[deprecated(
    "Acts::Experimental::PadBlueprintNode moved to Acts::PadBlueprintNode")]] =
    Acts::PadBlueprintNode;
}  // namespace Experimental

}  // namespace Acts
