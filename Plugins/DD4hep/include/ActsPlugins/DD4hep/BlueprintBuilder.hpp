// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintBuilder.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/Root/TGeoAxes.hpp"

#include <format>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <DD4hep/Detector.h>

namespace dd4hep {
class DetElement;
}  // namespace dd4hep

namespace Acts {
namespace Experimental {
class CylinderContainerBlueprintNode;
class CuboidContainerBlueprintNode;
}  // namespace Experimental
}  // namespace Acts

namespace ActsPlugins::DD4hep {

class PlaneLayerHelper;

// @TODO: Flip this
using GeometryAxes = ActsPlugins::TGeoAxes;

class DD4hepBackend {
 public:
  using Element = dd4hep::DetElement;
  using AxisDefinition = TGeoAxes;
  struct LayerSpec {
    std::optional<AxisDefinition> axes;
    std::optional<AxisDefinition> layerAxes;
    std::optional<std::string> layerName;
  };
  using DetectorElement = DD4hepDetectorElement;
  using DetectorElementPtr = std::shared_ptr<DetectorElement>;
  using DetectorElementFactory = std::function<DetectorElementPtr(
      const Element& detElement, AxisDefinition axes, double lengthScale)>;

  static DetectorElementPtr defaultElementFactory(const Element& detElement,
                                                  AxisDefinition axes,
                                                  double lengthScale);

  struct Config {
    DetectorElementFactory elementFactory = defaultElementFactory;
    const dd4hep::Detector* dd4hepDetector;
    double lengthScale = 1.0;
    std::reference_wrapper<const Acts::GeometryContext> gctx;
  };

  explicit DD4hepBackend(const Config& cfg, const Acts::Logger& logger);

  DetectorElementPtr createDetectorElement(const Element& detElement,
                                           AxisDefinition axes) const;

  std::vector<std::shared_ptr<Acts::Surface>> makeSurfaces(
      std::span<const Element> sensitives, const LayerSpec& layerSpec) const;

  std::optional<Acts::Transform3> lookupLayerTransform(
      const Element& element, const LayerSpec& layerSpec) const;

  std::shared_ptr<Acts::Experimental::StaticBlueprintNode> makeBeampipe() const;

  Element world() const;
  std::string nameOf(const Element& element) const;
  std::vector<Element> children(const Element& parent) const;
  Element parent(const Element& element) const;

  bool isSensitive(const Element& element) const;
  bool isBarrel(const Element& element) const;
  bool isEndcap(const Element& element) const;
  bool isTracker(const Element& element) const;

  /// Retrieves a named integer constant from the DD4hep detector description.
  /// The name is constructed by formatting @p fmt with @p args.
  template <typename... Args>
  int constant(std::format_string<Args...> fmt, Args&&... args) const {
    return m_cfg.dd4hepDetector->constant<int>(
        std::format(fmt, std::forward<Args>(args)...));
  }

  const Acts::Logger& logger() const { return *m_logger; }

 private:
  Config m_cfg;
  const Acts::Logger* m_logger;
};

using BlueprintBuilder = Acts::Experimental::BlueprintBuilder<DD4hepBackend>;
using ElementLayerAssembler =
    Acts::Experimental::ElementLayerAssembler<DD4hepBackend>;
using SensorLayerAssembler =
    Acts::Experimental::SensorLayerAssembler<DD4hepBackend>;
using SensorLayer = Acts::Experimental::SensorLayer<DD4hepBackend>;
using BarrelEndcapAssembler =
    Acts::Experimental::BarrelEndcapAssembler<DD4hepBackend>;

// plane
class PlaneLayerHelper {
 public:
  using Customizer = LayerHelper::Customizer;

  explicit PlaneLayerHelper(const BlueprintBuilder& builder)
      : m_builder{&builder} {}

  PlaneLayerHelper& setAxes(const std::string& axes) {
    m_axes = axes;
    return *this;
  }

  PlaneLayerHelper& setLayerAxes(const std::string& layerAxes) {
    m_layerAxes = layerAxes;
    return *this;
  }

  PlaneLayerHelper& setPattern(const std::string& pattern) {
    return setPattern(std::regex{pattern});
  }

  PlaneLayerHelper& setPattern(const std::regex& pattern) {
    m_pattern = pattern;
    return *this;
  }

  PlaneLayerHelper& setContainer(const dd4hep::DetElement& container) {
    m_container = container;
    return *this;
  }

  PlaneLayerHelper& setContainer(const std::string& name);

  PlaneLayerHelper& setEnvelope(const Acts::ExtentEnvelope& envelope) {
    m_envelope = envelope;
    return *this;
  }

  PlaneLayerHelper& setEmptyOk(bool emptyOk) {
    m_emptyOk = emptyOk;
    return *this;
  }

  PlaneLayerHelper& customize(Customizer customizer) {
    m_customizer = std::move(customizer);
    return *this;
  }

  auto& setAttachmentStrategy(
      std::optional<Acts::VolumeAttachmentStrategy> strategy) {
    m_attachmentStrategy = strategy;
    return *this;
  }

  // cuboid node
  std::shared_ptr<Acts::Experimental::CuboidContainerBlueprintNode> build()
      const;

  void addTo(Acts::Experimental::BlueprintNode& node) const;

 private:
  const BlueprintBuilder* m_builder;
  std::optional<std::string> m_axes;
  std::optional<std::string> m_layerAxes;
  std::optional<std::regex> m_pattern;
  std::optional<dd4hep::DetElement> m_container;
  std::optional<Acts::ExtentEnvelope> m_envelope;
  std::optional<Acts::VolumeAttachmentStrategy> m_attachmentStrategy;
  bool m_emptyOk = false;

  Customizer m_customizer;
};

}  // namespace ActsPlugins::DD4hep

// Explicit instantiation: suppress implicit instantiation in TUs that include
// this header. Definitions are instantiated in BlueprintBuilder.cpp.
// Placed at global scope so we open ::Acts::Experimental, not a nested Acts.
namespace Acts::Experimental {
extern template class BlueprintBuilder<ActsPlugins::DD4hep::DD4hepBackend>;
extern template class ElementLayerAssembler<ActsPlugins::DD4hep::DD4hepBackend>;
extern template class SensorLayerAssembler<ActsPlugins::DD4hep::DD4hepBackend>;
extern template class SensorLayer<ActsPlugins::DD4hep::DD4hepBackend>;
extern template class BarrelEndcapAssembler<ActsPlugins::DD4hep::DD4hepBackend>;
}  // namespace Acts::Experimental
