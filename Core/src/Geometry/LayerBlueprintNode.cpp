// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/LayerBlueprintNode.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/GraphViz.hpp"

namespace Acts::Experimental {

namespace detail {
struct LayerBlueprintNodeImpl {
  using LayerType = LayerBlueprintNode::LayerType;

  std::string m_name;

  std::vector<std::shared_ptr<Surface>> m_surfaces{};

  /// If a proto layer is already given externally, this node will not perform
  /// sizing from surfaces
  std::optional<MutableProtoLayer> m_protoLayer;

  Transform3 m_transform = Transform3::Identity();
  ExtentEnvelope m_envelope = ExtentEnvelope::Zero();
  LayerType m_layerType = LayerType::Cylinder;
  std::array<bool, 3> m_useCenterOfGravity = {true, true, true};
};
}  // namespace detail

LayerBlueprintNode::LayerBlueprintNode(std::string_view name)
    : StaticBlueprintNode(nullptr) {
  m_impl = std::make_unique<detail::LayerBlueprintNodeImpl>();
  m_impl->m_name = name;
}

LayerBlueprintNode::~LayerBlueprintNode() = default;

detail::LayerBlueprintNodeImpl& LayerBlueprintNode::impl() {
  assert(m_impl != nullptr);
  return *m_impl;
}

const detail::LayerBlueprintNodeImpl& LayerBlueprintNode::impl() const {
  assert(m_impl != nullptr);
  return *m_impl;
}

Volume& LayerBlueprintNode::build(const BlueprintOptions& options,
                                  const GeometryContext& gctx,
                                  const Logger& logger) {
  if (impl().m_surfaces.empty()) {
    ACTS_ERROR("LayerBlueprintNode: no surfaces provided");
    throw std::invalid_argument("LayerBlueprintNode: no surfaces provided");
  }

  ACTS_DEBUG(prefix() << "Building Layer " << name() << " from "
                      << impl().m_surfaces.size() << " surfaces");
  ACTS_VERBOSE(prefix() << " -> layer type: " << impl().m_layerType);
  ACTS_VERBOSE(prefix() << " -> transform:\n" << impl().m_transform.matrix());

  Extent extent;

  if (!impl().m_protoLayer.has_value()) {
    impl().m_protoLayer.emplace(gctx, impl().m_surfaces,
                                impl().m_transform.inverse());
    ACTS_VERBOSE(prefix() << "Built proto layer: "
                          << impl().m_protoLayer.value());
  } else {
    ACTS_VERBOSE(prefix() << "Using provided proto layer");
  }

  auto& protoLayer = impl().m_protoLayer.value();
  extent.addConstrain(protoLayer.extent, impl().m_envelope);

  ACTS_VERBOSE(prefix() << " -> layer extent: " << extent);

  buildVolume(extent, logger);
  assert(m_volume != nullptr && "Volume not built from proto layer");

  for (auto& surface : impl().m_surfaces) {
    m_volume->addSurface(surface);
  }

  return StaticBlueprintNode::build(options, gctx, logger);
}

void LayerBlueprintNode::buildVolume(const Extent& extent,
                                     const Logger& logger) {
  ACTS_VERBOSE(prefix() << "Building volume for layer " << name());
  using enum AxisDirection;
  using enum LayerType;

  std::shared_ptr<VolumeBounds> bounds;
  switch (impl().m_layerType) {
    case Cylinder:
    case Disc: {
      double minR = extent.min(AxisR);
      double maxR = extent.max(AxisR);
      double hlZ = extent.interval(AxisZ) / 2.0;
      bounds = std::make_shared<CylinderVolumeBounds>(minR, maxR, hlZ);
      break;
    }
    case Plane: {
      double hlX = extent.interval(AxisX) / 2.0;
      double hlY = extent.interval(AxisY) / 2.0;
      double hlZ = extent.interval(AxisZ) / 2.0;
      bounds = std::make_shared<CuboidVolumeBounds>(hlX, hlY, hlZ);
      break;
    }
  }

  assert(bounds != nullptr);

  ACTS_VERBOSE(prefix() << " -> bounds: " << *bounds);

  Transform3 transform = impl().m_transform;
  Vector3 translation = Vector3::Zero();
  if (impl().m_useCenterOfGravity.at(toUnderlying(AxisX))) {
    translation.x() = extent.medium(AxisX);
  }
  if (impl().m_useCenterOfGravity.at(toUnderlying(AxisY))) {
    translation.y() = extent.medium(AxisY);
  }
  if (impl().m_useCenterOfGravity.at(toUnderlying(AxisZ))) {
    translation.z() = extent.medium(AxisZ);
  }

  transform.translation() = translation;

  ACTS_VERBOSE(prefix() << " -> adjusted transform:\n" << transform.matrix());

  m_volume = std::make_unique<TrackingVolume>(transform, std::move(bounds),
                                              impl().m_name);
}

const std::string& LayerBlueprintNode::name() const {
  return impl().m_name;
}

LayerBlueprintNode& LayerBlueprintNode::setSurfaces(
    std::vector<std::shared_ptr<Surface>> surfaces) {
  impl().m_surfaces = std::move(surfaces);
  impl().m_protoLayer.reset();
  return *this;
}

const std::vector<std::shared_ptr<Surface>>& LayerBlueprintNode::surfaces()
    const {
  return impl().m_surfaces;
}

LayerBlueprintNode& LayerBlueprintNode::setProtoLayer(
    std::optional<MutableProtoLayer> protoLayer) {
  impl().m_protoLayer = std::move(protoLayer);
  impl().m_surfaces.clear();
  // also take ownership of the surfaces now
  for (auto& surface : impl().m_protoLayer.value().surfaces()) {
    impl().m_surfaces.push_back(surface->getSharedPtr());
  }
  return *this;
}

const MutableProtoLayer* LayerBlueprintNode::protoLayer() const {
  return impl().m_protoLayer.has_value() ? &impl().m_protoLayer.value()
                                         : nullptr;
}

LayerBlueprintNode& LayerBlueprintNode::setTransform(
    const Transform3& transform) {
  impl().m_transform = transform;
  return *this;
}

const Transform3& LayerBlueprintNode::transform() const {
  return impl().m_transform;
}

LayerBlueprintNode& LayerBlueprintNode::setEnvelope(
    const ExtentEnvelope& envelope) {
  impl().m_envelope = envelope;
  return *this;
}

const ExtentEnvelope& LayerBlueprintNode::envelope() const {
  return impl().m_envelope;
}

LayerBlueprintNode& LayerBlueprintNode::setLayerType(LayerType layerType) {
  impl().m_layerType = layerType;
  return *this;
}

const LayerBlueprintNode::LayerType& LayerBlueprintNode::layerType() const {
  return impl().m_layerType;
}

LayerBlueprintNode& LayerBlueprintNode::setUseCenterOfGravity(bool x, bool y,
                                                              bool z) {
  impl().m_useCenterOfGravity = {x, y, z};
  return *this;
}

void LayerBlueprintNode::addToGraphviz(std::ostream& os) const {
  std::stringstream ss;
  ss << "<br/><b>" + name() + "</b>";
  ss << "<br/>Layer";
  ss << "<br/><i>" << impl().m_layerType << "</i>";

  GraphViz::Node node{
      .id = name(), .label = ss.str(), .shape = GraphViz::Shape::Diamond};

  os << node;

  BlueprintNode::addToGraphviz(os);
}

}  // namespace Acts::Experimental
