// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TransformRange.hpp"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace Acts {

class Volume;
class TrackingVolume;
class VolumeBounds;
class PortalShellBase;
class CylinderContainerBlueprintNode;
class MaterialDesignatorBlueprintNode;
class StaticBlueprintNode;
class LayerBlueprintNode;

/// Base class for all nodes in the blueprint tree. This class defines the
/// three-phase construction process. The three phases are
///
/// -# **Build**: Construct volume representation + compute final sizing
/// -# **Connect**: Create and connect portals at volume boundaries
/// -# **Finalize**: Register portals with volumes + create acceleration
///    structures
///
/// During the *build* phase, the `build` method of all nodes in the tree are
/// called recursively. Some nodes, like @c CylinderContainerBlueprintNode,
/// will take action on the volumes returns from its children, and perform
/// sizing to connect them. See the @c CylinderContainerBlueprintNode and @c
/// CylinderVolumeStack documentation for details on how the sizing is carried
/// out.
class BlueprintNode {
 public:
  /// Can be default constructed
  BlueprintNode() = default;

  /// Virtual destructor to ensure correct cleanup
  virtual ~BlueprintNode() = default;

  /// Get the name of this node
  virtual const std::string& name() const = 0;

  /// This method is called during the *build* phase of the blueprint tree
  /// construction. It returns a single @c Volume which represents transform
  /// and bounds of the entire subtree. This does not have to correspond to the
  /// final @c TrackingVolume, some node types will produce temporary volume
  /// representations. Lifetime of the returned volume is managed by the source
  /// node!
  /// Nodes higher in the hierarchy will issue resizes down the tree hierarchy.
  /// This is not done through a direct hierarchy, but coordinated by the
  /// respective node type, by internally consulting its children.
  ///
  /// @note Generally, you should not need to to call this method directly.
  ///       The construction should usually be done through the special
  ///       @c Blueprint class.
  ///
  /// @param options The global construction options
  /// @param gctx The geometry context for construction (usually nominal)
  /// @param logger The logger to use for output during construction
  /// @return The volume used for communicating transform and size up the hierarchy
  virtual Volume& build(const BlueprintOptions& options,
                        const GeometryContext& gctx,
                        const Logger& logger = Acts::getDummyLogger()) = 0;

  /// This method is called during the *connect* phase. This phase handles the
  /// creation and connection of *portals* (instances of @c PortalLinkBase).
  /// After the build-phase has completed, the volume sizes are **final**. Each
  /// node will consult its fully sized volume to produce *boundary surfaces*.
  /// Each boundary surface is then turned into a @c TrivialPortalLink, which
  /// in turn produces a one-sided portal (see @c Portal documentation)
  ///
  /// Some nodes (like @c CylinderContainerBlueprintNode) will take action on
  /// their children, and unify the connected portals.
  ///
  /// After a node's processing has completed, it returns a reference to a @c
  /// PortalShellBase, which represents a set of portals in a specific geometry
  /// arrangement. The returned object lifetime is managed by the returning
  /// node.
  ///
  /// @param options The global construction options
  /// @param gctx The geometry context for construction (usually nominal)
  /// @param logger The logger to use for output during construction
  virtual PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) = 0;

  /// This method is called during the *finalize* phase. This phase handles:
  ///
  /// - Registering portals into their final volumes
  /// - Registering volumes into their parents
  /// - Creating navigation policies
  /// - (In future) Handle geometry identification assignment
  ///
  /// At the end of this phase, each node will have transfered any temporary
  /// resources created during the build, that need to be retained, into the
  /// final @c TrackingGeometry, and can be safely destroyed.
  ///
  /// @note The parent for volumes, portals, etc to be registered in is passed in as an
  ///       argument, rather than being implicitly determined from the **parent
  ///       node**. This is done so that nodes can remove themselves from the
  ///       final volume hierarchy, like container nodes or the
  ///       @c MaterialDesignatorBlueprintNode.
  ///
  /// @param options The global construction options
  /// @param gctx The geometry context for construction (usually nominal)
  /// @param parent The parent volume to register in
  /// @param logger The logger to use for output during construction
  virtual void finalize(const BlueprintOptions& options,
                        const GeometryContext& gctx, TrackingVolume& parent,
                        const Logger& logger = Acts::getDummyLogger()) = 0;

  StaticBlueprintNode& addStaticVolume(
      std::unique_ptr<TrackingVolume> volume,
      const std::function<void(StaticBlueprintNode& cylinder)>& callback = {});

  StaticBlueprintNode& addStaticVolume(
      const Transform3& transform, std::shared_ptr<VolumeBounds> volbounds,
      const std::string& volumeName = "undefined",
      const std::function<void(StaticBlueprintNode& cylinder)>& callback = {});

  CylinderContainerBlueprintNode& addCylinderContainer(
      const std::string& name, BinningValue direction,
      const std::function<void(CylinderContainerBlueprintNode& cylinder)>&
          callback = {});

  MaterialDesignatorBlueprintNode& addMaterial(
      const std::string& name,
      const std::function<void(MaterialDesignatorBlueprintNode& material)>&
          callback = {});

  LayerBlueprintNode& addLayer(
      const std::string& name,
      const std::function<void(LayerBlueprintNode& layer)>& callback = {});

  BlueprintNode& addChild(std::shared_ptr<BlueprintNode> child);

  using MutableChildRange =
      detail::TransformRange<detail::Dereference,
                             std::vector<std::shared_ptr<BlueprintNode>>>;

  using ChildRange =
      detail::TransformRange<detail::ConstDereference,
                             const std::vector<std::shared_ptr<BlueprintNode>>>;

  MutableChildRange children();
  ChildRange children() const;

  void clearChildren();

  std::size_t depth() const;

  void graphViz(std::ostream& os) const;
  virtual void addToGraphviz(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const BlueprintNode& node) {
    node.toStream(os);
    return os;
  }

 protected:
  virtual void toStream(std::ostream& os) const;

  std::string prefix() const;
  std::string indent() const;

 private:
  std::size_t m_depth{0};
  std::vector<std::shared_ptr<BlueprintNode>> m_children{};
};
};  // namespace Acts
