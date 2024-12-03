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

class BlueprintNode {
 public:
  BlueprintNode() = default;

  virtual ~BlueprintNode() = default;

  virtual const std::string& name() const = 0;

  virtual Volume& build(const BlueprintOptions& options,
                        const GeometryContext& gctx,
                        const Logger& logger = Acts::getDummyLogger()) = 0;

  virtual PortalShellBase& connect(
      const BlueprintOptions& options, const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) = 0;

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
