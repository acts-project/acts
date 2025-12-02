// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/DD4hep/LayerBlueprintNode.hpp"

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

namespace ActsPlugins::DD4hep {

std::shared_ptr<DD4hepDetectorElement> BlueprintBuilder::defaultElementFactory(
    const dd4hep::DetElement& detElement, const std::string& axes,
    double lengthScale) {
  return std::make_shared<DD4hepDetectorElement>(detElement, axes, lengthScale);
}

std::shared_ptr<DD4hepDetectorElement> BlueprintBuilder::createDetectorElement(
    const dd4hep::DetElement& detElement, const std::string& axes) const {
  auto elem = m_cfg.elementFactory(detElement, axes, m_cfg.lengthScale);

  detElement.addExtension<DD4hepDetectorElementExtension>(
      new dd4hep::rec::StructExtension(DD4hepDetectorElementExtension(elem)));

  return elem;
}

namespace {
void visitSubtree(
    const dd4hep::DetElement& detElement,
    const std::function<void(const dd4hep::DetElement&)>& visitor) {
  visitor(detElement);

  for (auto& [name, child] : detElement.children()) {
    visitSubtree(child, visitor);
  }
}

}  // namespace

std::optional<dd4hep::DetElement> BlueprintBuilder::findDetElementByName(
    const dd4hep::DetElement& parent, const std::string& name) {
  if (parent.name() == name) {
    return parent;
  }

  for (const auto& [childName, child] : parent.children()) {
    auto result = findDetElementByName(child, name);
    if (result.has_value()) {
      return result;
    }
  }

  return std::nullopt;
}

std::optional<dd4hep::DetElement> BlueprintBuilder::findDetElementByName(
    const std::string& name) {
  return findDetElementByName(world(), name);
}

std::vector<dd4hep::DetElement> BlueprintBuilder::findDetElementByNamePattern(
    const dd4hep::DetElement& parent, const std::regex& pattern) {
  std::vector<dd4hep::DetElement> matches;

  visitSubtree(parent, [&](const auto& elem) {
    if (std::regex_match(elem.name(), pattern)) {
      matches.push_back(elem);
    }
  });

  return matches;
}

dd4hep::DetElement BlueprintBuilder::world() const {
  return m_cfg.dd4hepDetector->world();
}

std::vector<dd4hep::DetElement> BlueprintBuilder::resolveSensitives(
    const dd4hep::DetElement& detElement) {
  std::vector<dd4hep::DetElement> sensitives;
  visitSubtree(detElement, [&](const auto& elem) {
    if (elem.volume().isSensitive()) {
      sensitives.push_back(elem);
    }
  });
  return sensitives;
}

std::shared_ptr<LayerBlueprintNode> BlueprintBuilder::addLayer(
    const dd4hep::DetElement& detElement, const std::string& axes) {
  auto node = std::make_shared<LayerBlueprintNode>(detElement.name());

  auto sensitives = resolveSensitives(detElement);

  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  surfaces.reserve(sensitives.size());

  for (const auto& sensitive : sensitives) {
    auto elem = createDetectorElement(sensitive, axes);
    surfaces.push_back(elem->surface().getSharedPtr());
  }

  node->setSurfaces(std::move(surfaces));

  // @TODO: Try to auto-detect what kind of layer this is

  return node;
}

std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
BlueprintBuilder::addLayers(const dd4hep::DetElement& container,
                            const std::string& axes,
                            Acts::AxisDirection direction,
                            const std::regex& layerPattern) {
  using enum Acts::AxisDirection;
  using namespace Acts::UnitLiterals;
  auto node =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          container.name(), direction);
  for (const auto& element :
       findDetElementByNamePattern(container, layerPattern)) {
    auto layer = addLayer(element, axes);
    layer->setEnvelope(Acts::ExtentEnvelope{}
                           .set(AxisZ, {-5_mm, 5_mm})
                           .set(AxisR, {-5_mm, 5_mm}));

    node->addChild(layer);
  }

  return node;
}

}  // namespace ActsPlugins::DD4hep
