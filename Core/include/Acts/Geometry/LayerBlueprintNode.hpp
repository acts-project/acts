// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/StaticBlueprintNode.hpp"

#include <iosfwd>

namespace Acts {

class LayerBlueprintNode : public StaticBlueprintNode {
 public:
  enum class LayerType { Cylinder, Disc, Plane };

  LayerBlueprintNode(const std::string& name)
      : StaticBlueprintNode{nullptr}, m_name(name) {}

  const std::string& name() const override;

  Volume& build(const Options& options, const GeometryContext& gctx,
                const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(const Options& options, const GeometryContext& gctx,
                TrackingVolume& parent,
                const Logger& logger = Acts::getDummyLogger()) override;

  LayerBlueprintNode& setSurfaces(
      std::vector<std::shared_ptr<Surface>> surfaces);

  const std::vector<std::shared_ptr<Surface>>& surfaces() const;

  LayerBlueprintNode& setTransform(const Transform3& transform);

  const Transform3& transform() const;

  LayerBlueprintNode& setEnvelope(const ExtentEnvelope& envelope);

  const ExtentEnvelope& envelope() const;

  LayerBlueprintNode& setLayerType(LayerType layerType);

  const LayerType& layerType() const;

  LayerBlueprintNode& setNavigationPolicyFactory(
      std::shared_ptr<NavigationPolicyFactory> navigationPolicyFactory)
      override;

  const NavigationPolicyFactory* navigationPolicyFactory() const;

  void addToGraphviz(std::ostream& os) const override;

 private:
  void buildVolume(const Extent& extent, const Logger& logger);

  std::string m_name;
  std::vector<std::shared_ptr<Surface>> m_surfaces{};
  Transform3 m_transform = Transform3::Identity();
  ExtentEnvelope m_envelope = ExtentEnvelope::Zero();
  LayerType m_layerType = LayerType::Cylinder;
};

std::ostream& operator<<(std::ostream& os, LayerBlueprintNode::LayerType type);

}  // namespace Acts
