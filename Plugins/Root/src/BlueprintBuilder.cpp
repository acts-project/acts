// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/BlueprintBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/detail/BlueprintBuilder_impl.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/Root/TGeoSurfaceConverter.hpp"

#include <cstddef>
#include <functional>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"

namespace ActsPlugins {

TGeoBlueprintBuilderBackend::DetectorElementPtr
TGeoBlueprintBuilderBackend::defaultElementFactory(
    const TGeoDetectorElement::Identifier& identifier, const TGeoNode& tGeoNode,
    const TGeoMatrix& tGeoMatrix, AxisDefinition axes, double lengthScale,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  return std::make_shared<TGeoDetectorElement>(
      identifier, tGeoNode, tGeoMatrix, axes, lengthScale, std::move(material));
}

TGeoBlueprintBuilderBackend::TGeoBlueprintBuilderBackend(
    const Config& cfg, const Acts::Logger& logger)
    : m_cfg(cfg), m_logger(&logger) {
  if (m_cfg.root == nullptr) {
    throw std::invalid_argument(
        "TGeoBlueprintBuilderBackend: root node is null");
  }
  m_world = makeElement(*m_cfg.root, nullptr);
}

auto TGeoBlueprintBuilderBackend::makeElement(
    const TGeoNode& node, std::shared_ptr<const NodeContext> parent) const
    -> Element {
  return Element{std::make_shared<NodeContext>(
      NodeContext{.node = &node, .parent = std::move(parent)})};
}

const TGeoBlueprintBuilderBackend::NodeContext&
TGeoBlueprintBuilderBackend::contextOf(const Element& element) const {
  if (element.context == nullptr) {
    throw std::invalid_argument(
        "TGeoBlueprintBuilderBackend: invalid element handle");
  }
  return *element.context;
}

TGeoDetectorElement::Identifier
TGeoBlueprintBuilderBackend::defaultIdentifierFor(
    const Element& element) const {
  const auto path = pathOf(element);
  return static_cast<TGeoDetectorElement::Identifier>(
      std::hash<std::string>{}(path));
}

TGeoBlueprintBuilderBackend::DetectorElementPtr
TGeoBlueprintBuilderBackend::createDetectorElement(const Element& element,
                                                   AxisDefinition axes) const {
  const auto& context = contextOf(element);
  const auto identifier = m_cfg.identifierProvider != nullptr
                              ? m_cfg.identifierProvider(element)
                              : defaultIdentifierFor(element);
  const auto transform = transformOf(element);

  auto detectorElement = m_cfg.elementFactory(
      identifier, *context.node, transform, axes, m_cfg.lengthScale, nullptr);
  m_detectorElementStore.push_back(detectorElement);
  return detectorElement;
}

std::vector<std::shared_ptr<Acts::Surface>>
TGeoBlueprintBuilderBackend::makeSurfaces(std::span<const Element> sensitives,
                                          const LayerSpec& layerSpec) const {
  if (!layerSpec.axes.has_value()) {
    throw std::runtime_error(
        "TGeoBlueprintBuilderBackend::makeSurfaces: axes not set");
  }

  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  surfaces.reserve(sensitives.size());

  for (const auto& sensitive : sensitives) {
    auto detectorElement =
        createDetectorElement(sensitive, layerSpec.axes.value());
    surfaces.push_back(detectorElement->surface().getSharedPtr());
  }

  return surfaces;
}

std::optional<Acts::Transform3>
TGeoBlueprintBuilderBackend::lookupLayerTransform(
    const Element& element, const LayerSpec& layerSpec) const {
  if (!layerSpec.layerAxes.has_value()) {
    return std::nullopt;
  }

  const auto& context = contextOf(element);
  return TGeoSurfaceConverter::transformFromShape(
      *context.node->GetVolume()->GetShape(), transformOf(element),
      layerSpec.layerAxes.value(), m_cfg.lengthScale);
}

TGeoBlueprintBuilderBackend::Element TGeoBlueprintBuilderBackend::world()
    const {
  return m_world;
}

std::string TGeoBlueprintBuilderBackend::nameOf(const Element& element) const {
  const auto& context = contextOf(element);
  if (const auto* volume = context.node->GetVolume();
      volume != nullptr && volume->GetName() != nullptr) {
    return volume->GetName();
  }
  if (context.node->GetName() != nullptr) {
    return context.node->GetName();
  }
  return {};
}

std::vector<TGeoBlueprintBuilderBackend::Element>
TGeoBlueprintBuilderBackend::children(const Element& parent) const {
  const auto& context = contextOf(parent);
  std::vector<Element> result;
  result.reserve(context.node->GetNdaughters());
  for (int i = 0; i < context.node->GetNdaughters(); ++i) {
    const TGeoNode* child = context.node->GetDaughter(i);
    if (child == nullptr) {
      continue;
    }
    result.push_back(makeElement(*child, parent.context));
  }
  return result;
}

TGeoBlueprintBuilderBackend::Element TGeoBlueprintBuilderBackend::parent(
    const Element& element) const {
  return Element{contextOf(element).parent};
}

bool TGeoBlueprintBuilderBackend::isSensitive(const Element& element) const {
  return m_cfg.sensitivePredicate != nullptr &&
         m_cfg.sensitivePredicate(element);
}

const TGeoNode& TGeoBlueprintBuilderBackend::nodeOf(
    const Element& element) const {
  return *contextOf(element).node;
}

TGeoHMatrix TGeoBlueprintBuilderBackend::transformOf(
    const Element& element) const {
  std::vector<const NodeContext*> chain;
  for (const NodeContext* current = element.context.get(); current != nullptr;
       current = current->parent.get()) {
    chain.push_back(current);
  }

  std::optional<TGeoHMatrix> transform = std::nullopt;
  for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
    const TGeoMatrix* localMatrix = (*it)->node->GetMatrix();
    if (localMatrix == nullptr) {
      continue;
    }
    if (!transform.has_value()) {
      transform.emplace(*localMatrix);
      continue;
    }
    transform = TGeoCombiTrans(*transform) * TGeoCombiTrans(*localMatrix);
  }
  return transform.value_or(TGeoHMatrix{});
}

std::string TGeoBlueprintBuilderBackend::pathOf(const Element& element) const {
  std::vector<std::string> names;
  for (const NodeContext* current = element.context.get(); current != nullptr;
       current = current->parent.get()) {
    if (current->node->GetName() != nullptr) {
      names.emplace_back(current->node->GetName());
    } else {
      names.emplace_back(current->node->GetVolume()->GetName());
    }
  }

  if (names.empty()) {
    return {};
  }

  std::string path = names.back();
  for (auto it = std::next(names.rbegin()); it != names.rend(); ++it) {
    path += "|";
    path += *it;
  }
  return path;
}

}  // namespace ActsPlugins

namespace Acts::Experimental {
template class BlueprintBuilder<ActsPlugins::TGeoBlueprintBuilderBackend>;
template class ElementLayerAssembler<ActsPlugins::TGeoBlueprintBuilderBackend>;
template class SensorLayerAssembler<ActsPlugins::TGeoBlueprintBuilderBackend>;
template class SensorLayer<ActsPlugins::TGeoBlueprintBuilderBackend>;
}  // namespace Acts::Experimental
