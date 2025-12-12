// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
// Needed for explicit instantiation of template methods.
#include "Acts/Geometry/detail/BlueprintBuilder_impl.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/Root/TGeoSurfaceConverter.hpp"

#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <DD4hep/DetElement.h>
#include <DD4hep/DetType.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

namespace ActsPlugins::DD4hep {

using ActsPlugins::TGeoAxes;

DD4hepBackend::DetectorElementPtr DD4hepBackend::defaultElementFactory(
    const Element& detElement, AxisDefinition axes, double lengthScale) {
  return std::make_shared<DD4hepDetectorElement>(detElement, axes, lengthScale);
}

DD4hepBackend::DD4hepBackend(const Config& cfg, const Acts::Logger& logger)
    : m_cfg(cfg), m_logger(&logger) {
  if (m_cfg.dd4hepDetector == nullptr) {
    throw std::invalid_argument("DD4hepBackend: dd4hepDetector is null");
  }
}

DD4hepBackend::DetectorElementPtr DD4hepBackend::createDetectorElement(
    const Element& detElement, AxisDefinition axes) const {
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
    (void)name;
    visitSubtree(child, visitor);
  }
}

}  // namespace

DD4hepBackend::Element DD4hepBackend::world() const {
  return m_cfg.dd4hepDetector->world();
}

std::string DD4hepBackend::nameOf(const Element& element) const {
  return element.name();
}

std::vector<DD4hepBackend::Element> DD4hepBackend::children(
    const Element& parent) const {
  std::vector<Element> result;
  result.reserve(parent.children().size());
  for (const auto& [name, child] : parent.children()) {
    (void)name;
    result.push_back(child);
  }
  return result;
}

DD4hepBackend::Element DD4hepBackend::parent(const Element& element) const {
  return element.parent();
}

bool DD4hepBackend::isSensitive(const Element& element) const {
  return element.volume().isSensitive();
}

bool DD4hepBackend::isBarrel(const Element& element) const {
  return dd4hep::DetType{element.typeFlag()}.is(dd4hep::DetType::BARREL);
}

bool DD4hepBackend::isEndcap(const Element& element) const {
  return dd4hep::DetType{element.typeFlag()}.is(dd4hep::DetType::ENDCAP);
}

bool DD4hepBackend::isTracker(const Element& element) const {
  return dd4hep::DetType{element.typeFlag()}.is(dd4hep::DetType::TRACKER);
}

std::vector<std::shared_ptr<Acts::Surface>> DD4hepBackend::makeSurfaces(
    std::span<const dd4hep::DetElement> sensitives,
    const LayerSpec& layerSpec) const {
  if (!layerSpec.axes.has_value()) {
    throw std::runtime_error("DD4hepBackend::makeSurfaces: axes not set");
  }

  ACTS_DEBUG("Using " << sensitives.size() << " sensitive elements.");
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  surfaces.reserve(sensitives.size());

  for (const auto& sensitive : sensitives) {
    auto elem = createDetectorElement(sensitive, layerSpec.axes.value());
    surfaces.push_back(elem->surface().getSharedPtr());
  }

  return surfaces;
}

std::optional<Acts::Transform3> DD4hepBackend::lookupLayerTransform(
    const dd4hep::DetElement& element, const LayerSpec& layerSpec) const {
  if (layerSpec.layerAxes.has_value()) {
    ACTS_DEBUG("Finding layer transform automatically using layer axes: "
               << layerSpec.layerAxes.value());
    Acts::Transform3 layerTransform = TGeoSurfaceConverter::transformFromShape(
        *element.placement().ptr()->GetVolume()->GetShape(),
        element.nominal().worldTransformation(), layerSpec.layerAxes.value(),
        m_cfg.lengthScale);

    ACTS_VERBOSE(" -> Layer transform:\n" << layerTransform.matrix());
    return layerTransform;
  }

  return std::nullopt;
}

std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
DD4hepBackend::makeBeampipe() const {
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

  ACTS_INFO("-> Created beampipe volume: "
            << volume->volumeBounds() << " transform:\n"
            << volume->localToGlobalTransform(m_cfg.gctx.get()).matrix());

  return std::make_shared<Acts::Experimental::StaticBlueprintNode>(
      std::move(volume));
}

// container for plane
PlaneLayerHelper& PlaneLayerHelper::setContainer(const std::string& name) {
  m_container = m_builder->findDetElementByName(name);
  if (!m_container.has_value()) {
    throw std::runtime_error("Could not find DetElement with name " + name +
                             " in PlaneLayerHelper");
  }
  return *this;
}

// plane builder
std::shared_ptr<Acts::Experimental::CuboidContainerBlueprintNode>
PlaneLayerHelper::build() const {
  const auto& logger = m_builder->logger();

  if (!m_axes.has_value()) {
    throw std::runtime_error("Axes not set in PlaneLayerHelper");
  }

  if (!m_pattern.has_value()) {
    throw std::runtime_error("Pattern not set in PlaneLayerHelper");
  }

  if (!m_container.has_value()) {
    throw std::runtime_error("Container not set in PlaneLayerHelper");
  }

  const auto& container = m_container.value();
  ACTS_INFO("[PlaneHelper]: building plane container " << container.name());

  auto node =
      std::make_shared<Acts::Experimental::CuboidContainerBlueprintNode>(
          container.name(), Acts::AxisDirection::AxisZ);

  if (m_attachmentStrategy.has_value()) {
    node->setAttachmentStrategy(m_attachmentStrategy.value());
  }

  const auto& layerElements =
      m_builder->findDetElementByNamePattern(container, m_pattern.value());

  if (layerElements.empty()) {
    ACTS_LOG(m_emptyOk ? Acts::Logging::INFO : Acts::Logging::ERROR,
             "No layers found in container " << container.name()
                                             << " matching pattern (plane)");
    if (!m_emptyOk) {
      throw std::runtime_error(
          std::format("No layers found in container {} matching pattern",
                      container.name()));
    }
  }

  for (const auto& element : layerElements) {
    auto layer =
        m_builder->makeLayer(element, m_axes.value(), m_layerAxes, m_layerType);
    if (m_envelope.has_value()) {
      layer->setEnvelope(m_envelope.value());
    }

    if (m_customizer) {
      layer = m_customizer(element, std::move(layer));
    }

    node->addChild(layer);
  }

  return node;
}  // plane builder

}  // namespace ActsPlugins::DD4hep

// Explicit template instantiation for DD4hepBackend. Ensures all template
// code is compiled in this TU; other TUs use extern template and link here.
// Must be in ::Acts::Experimental (at global scope) to match the template defs.
namespace Acts::Experimental {
template class BlueprintBuilder<ActsPlugins::DD4hep::DD4hepBackend>;
template class ElementLayerAssembler<ActsPlugins::DD4hep::DD4hepBackend>;
template class SensorLayerAssembler<ActsPlugins::DD4hep::DD4hepBackend>;
template class SensorLayer<ActsPlugins::DD4hep::DD4hepBackend>;
template class BarrelEndcapAssembler<ActsPlugins::DD4hep::DD4hepBackend>;
}  // namespace Acts::Experimental
