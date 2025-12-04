// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/Root/TGeoSurfaceConverter.hpp"
#include "ActsPlugins/Root/TGeoVolumeConverter.hpp"

#include <DD4hep/DetElement.h>
#include <DD4hep/DetType.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <boost/algorithm/string/join.hpp>

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

std::string BlueprintBuilder::getPathToElementName(
    const dd4hep::DetElement& elem) const {
  std::vector<std::string> names;
  names.emplace_back(elem.name());
  auto parent = elem.parent();
  while (parent != world()) {
    names.emplace_back(parent.name());
    parent = parent.parent();
  }
  std::ranges::reverse(names);
  return boost::algorithm::join(names, "|");
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

namespace {
Acts::Transform3 convertTGeoTransform(const TGeoShape& shape,
                                      const TGeoMatrix& transform,
                                      const std::string& axes,
                                      double lengthScale) {
  // This is somewhat duplicated from the TGeoDetectorElement constructor
  //
  const Double_t* translation = transform.GetTranslation();
  const Double_t* rotation = transform.GetRotationMatrix();

  auto [cBounds, cTransform, cThickness] =
      ActsPlugins::TGeoSurfaceConverter::cylinderComponents(
          shape, rotation, translation, axes, lengthScale);
  if (cBounds != nullptr) {
    return cTransform;
  }

  auto [dBounds, dTransform, dThickness] = TGeoSurfaceConverter::discComponents(
      shape, rotation, translation, axes, lengthScale);
  if (dBounds != nullptr) {
    return dTransform;
  }

  auto [pBounds, pTransform, pThickness] =
      TGeoSurfaceConverter::planeComponents(shape, rotation, translation, axes,
                                            lengthScale);
  if (pBounds != nullptr) {
    return pTransform;
  }

  throw std::runtime_error(
      "Could not extract transform from TGeoShape of type " +
      std::string(shape.ClassName()));
}

}  // namespace

std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
BlueprintBuilder::makeLayer(const dd4hep::DetElement& detElement,
                            const std::string& axes,
                            std::optional<std::string> layerAxes) const {
  ACTS_DEBUG("Adding layer from element: " << detElement.name());
  auto node = std::make_shared<Acts::Experimental::LayerBlueprintNode>(
      getPathToElementName(detElement));

  auto sensitives = resolveSensitives(detElement);
  ACTS_DEBUG("  Found " << sensitives.size() << " sensitive elements.");

  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  surfaces.reserve(sensitives.size());

  for (const auto& sensitive : sensitives) {
    auto elem = createDetectorElement(sensitive, axes);
    surfaces.push_back(elem->surface().getSharedPtr());
  }

  node->setSurfaces(std::move(surfaces));

  if (layerAxes.has_value()) {
    ACTS_DEBUG("Finding layer transform automatically using layer axes: "
               << layerAxes.value());
    Acts::Transform3 layerTransform = convertTGeoTransform(
        *detElement.placement().ptr()->GetVolume()->GetShape(),
        detElement.nominal().worldTransformation(), layerAxes.value(),
        m_cfg.lengthScale);

    ACTS_VERBOSE(" -> Layer transform:\n" << layerTransform.matrix());
    node->setTransform(layerTransform);
  }

  // @TODO: Try to auto-detect what kind of layer this is

  return node;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
BlueprintBuilder::makeBeampipe() const {
  std::optional<dd4hep::DetElement> beampipeElement = std::nullopt;

  visitSubtree(world(), [&](const dd4hep::DetElement& elem) {
    dd4hep::DetType subDetType{elem.typeFlag()};
    if (subDetType.is(dd4hep::DetType::BEAMPIPE)) {
      if (beampipeElement.has_value()) {
        ACTS_WARNING("Multiple beampipe elements found, using first: "
                     << beampipeElement->name()
                     << ", ignoring: " << elem.name());
        return;
      }
      beampipeElement = elem;
    }
  });

  if (!beampipeElement.has_value()) {
    ACTS_ERROR("No beampipe element found in DD4hep detector.");
    throw std::runtime_error("No beampipe element found in DD4hep detector.");
  }

  ACTS_INFO("Beampipe element found: " << beampipeElement->name());

  const auto& tgTransform = beampipeElement->nominal().worldTransformation();
  auto [bounds, transform, thickness] =
      ActsPlugins::TGeoSurfaceConverter::cylinderComponents(
          *beampipeElement->placement().ptr()->GetVolume()->GetShape(),
          tgTransform.GetRotationMatrix(), tgTransform.GetTranslation(), "XYZ",
          m_cfg.lengthScale);

  if (bounds == nullptr) {
    ACTS_ERROR("Beampipe element shape could not be converted to cylinder.");
    throw std::runtime_error(
        "Beampipe element shape could not be converted to cylinder.");
  }

  std::shared_ptr volBounds = std::make_shared<Acts::CylinderVolumeBounds>(
      0, bounds->get(Acts::CylinderBounds::eR),
      bounds->get(Acts::CylinderBounds::eHalfLengthZ));

  std::unique_ptr volume = std::make_unique<Acts::TrackingVolume>(
      transform, volBounds, beampipeElement->name());

  ACTS_INFO("-> Created beampipe volume: " << volume->volumeBounds()
                                           << " transform:\n"
                                           << volume->transform().matrix());

  return std::make_shared<Acts::Experimental::StaticBlueprintNode>(
      std::move(volume));
}

std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
BlueprintBuilder::addLayers(const dd4hep::DetElement& container,
                            const std::string& axes,
                            Acts::AxisDirection direction,
                            const std::regex& layerPattern,
                            const Acts::ExtentEnvelope& envelope) {
  ACTS_DEBUG("Adding layers to container " << container.name());
  using enum Acts::AxisDirection;
  using namespace Acts::UnitLiterals;
  auto node =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          container.name(), direction);
  const auto& layerElements =
      findDetElementByNamePattern(container, layerPattern);

  if (container.children().empty()) {
    ACTS_WARNING("Container " << container.name()
                              << " has no children, no layers added.");
  }
  for (const auto& element : layerElements) {
    auto layer = makeLayer(element, axes);
    layer->setEnvelope(envelope);

    node->addChild(layer);
  }

  return node;
}

LayerHelper BlueprintBuilder::layerHelper() {
  return LayerHelper(*this);
}

std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
LayerHelper::build() const {
  const auto& logger = m_builder->logger();

  if (!m_layerType.has_value()) {
    throw std::runtime_error("Layer type not set in LayerHelper");
  }

  if (!m_axes.has_value()) {
    throw std::runtime_error("Axes not set in LayerHelper");
  }

  if (!m_pattern.has_value()) {
    throw std::runtime_error("Pattern not set in LayerHelper");
  }

  if (!m_container.has_value()) {
    throw std::runtime_error("Container not set in LayerHelper");
  }

  Acts::AxisDirection axisDir = m_layerType == LayerType::Cylinder
                                    ? Acts::AxisDirection::AxisR
                                    : Acts::AxisDirection::AxisZ;

  const auto& container = m_container.value();
  auto node =
      std::make_shared<Acts::Experimental::CylinderContainerBlueprintNode>(
          container.name(), axisDir);

  const auto& layerElements =
      m_builder->findDetElementByNamePattern(container, m_pattern.value());

  if (layerElements.empty()) {
    ACTS_LOG(m_emptyOk ? Acts::Logging::INFO : Acts::Logging::ERROR,
             "No layers found in container " << container.name()
                                             << " matching pattern");
    if (!m_emptyOk) {
      throw std::runtime_error(
          std::format("No layers found in container {} matching pattern",
                      container.name()));
    }
  }
  for (const auto& element : layerElements) {
    auto layer = m_builder->makeLayer(element, m_axes.value(), m_layerAxes);
    if (m_envelope.has_value()) {
      layer->setEnvelope(m_envelope.value());
    }
    node->addChild(layer);

    if (m_customizer) {
      m_customizer(element, *layer);
    }
  }

  return node;
}

}  // namespace ActsPlugins::DD4hep
