// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/IVolumeLookup.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <optional>

namespace Acts {

class GeometryContext;

/// This class is the top-level entry point to build a tracking geometry using
/// the blueprint building mechanism. It forms the root of a tree of nodes where
/// each node performs a portion of the construction. This top-level class has
/// the main construction methods that execute the construction of the geometry.
///
/// ```
///            +---------------+  +-----------+
///            |               |  |           |
///            |     Root      |  |           v
///            |               |  |   +---------------+
///            +---------------+  |   |               |
///                    |          |   |    Child 1    |
///         +----------+          |   |               |
///         v          +----------+   +---------------+
/// +---------------+                         |
/// |               |          +--------------+
/// |    Child 2    |          v         +----------+
/// |               |     .---------.    |          |
/// +---------------+    /           \   |          v
///                     (  Proc node  )  |  +---------------+
///                      `.         ,'   |  |               |
///                        `-------'     |  |    Child 3    |
///                            |         |  |               |
///                            |         |  +---------------+
///                            +---------+
/// ```
///
/// The construction phases are documented in @c BlueprintNode, which is the
/// base class for all nodes in the tree.
/// @note This class inherits from @c BlueprintNode, but hides the main
///       blueprint construction phase overloads. The @c Blueprint class is
///       only ever intended to be the top-level node, and not anywhere else
///       in the tree.
class Blueprint : public BlueprintNode {
 public:
  /// Configuration for building a blueprint tracking geometry.
  struct Config {
    /// Determine how much envelope space to produce from the highest volume
    /// in the geometry hierarchy and the world volume.
    ExtentEnvelope envelope = ExtentEnvelope::Zero();
    /// Apply a bound deduplication on the world volume. It ensures
    /// that equivalent bounds are instantiated only once & recycled
    /// across the geometry components
    bool boundDeduplication{true};
    /// Optional global placement applied after sizing, before portal connection
    /// and navigation-policy construction. Transforms are left-multiplied by
    /// this transform, including the world. Shared objects are placed once.
    /// External alignment placements, pre-existing portals/navigation policies/
    /// volume lookups,
    /// and legacy layer/volume arrays or dense volumes are unsupported and
    /// rejected before any placement is applied. Custom nodes must implement
    /// @ref BlueprintNode::visitBuiltVolumes.
    std::optional<Transform3> globalTransform = std::nullopt;
  };

  /// Constructor from a config object
  /// @param config The configuration object
  explicit Blueprint(const Config& config);

  /// Construct the tracking geometry from the blueprint tree
  /// @param options The construction options, see @c BlueprintOptions
  /// @param gctx The geometry context for construction. In almost all cases,
  ///             this should be the *nominal* geometry context
  /// @param logger The logger to use for output during construction
  /// @return Unique pointer to the constructed tracking geometry
  std::unique_ptr<TrackingGeometry> construct(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger());

  /// Set the lookup factory for the world only.
  /// An empty factory leaves the volume's lookup unchanged. Descendant policies
  /// are installed first, so subtree indices can stop at explicit overrides.
  /// @param factory Factory returning a non-null lookup
  /// @return This node for chaining
  Blueprint& setVolumeLookupFactory(VolumeLookupFactory factory);

  /// Return the optional lookup factory.
  /// @return Configured factory, which may be empty
  const VolumeLookupFactory& volumeLookupFactory() const;

 protected:
  /// The name of the blueprint node, always "Root"
  /// @return The name
  const std::string& name() const override;

  /// @copydoc BlueprintNode::build
  Volume& build(const BlueprintOptions& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// @copydoc BlueprintNode::connect
  PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  /// @copydoc BlueprintNode::finalize
  void finalize(const BlueprintOptions& options, const GeometryContext& gctx,
                TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  /// @copydoc BlueprintNode::addToGraphviz
  void addToGraphviz(std::ostream& os) const override;

 private:
  VolumeLookupFactory m_volumeLookupFactory;

  Config m_cfg;
};

namespace Experimental {
/// @deprecated The blueprint geometry moved out of the `Acts::Experimental`
///             namespace. Use @ref Acts::Blueprint instead. This alias is kept
///             for backward compatibility and will be removed.
using Blueprint
    [[deprecated("Acts::Experimental::Blueprint moved to Acts::Blueprint")]] =
        Acts::Blueprint;
}  // namespace Experimental

}  // namespace Acts
