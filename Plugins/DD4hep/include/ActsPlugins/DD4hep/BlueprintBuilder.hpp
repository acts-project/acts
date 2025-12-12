// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"

#include <functional>
#include <memory>
#include <optional>
#include <regex>
#include <string>
#include <utility>

namespace dd4hep {
class DetElement;
}

namespace Acts {
namespace Experimental {
class CylinderContainerBlueprintNode;
}  // namespace Experimental
}  // namespace Acts

namespace ActsPlugins {

namespace DD4hep {

class LayerBlueprintNode;
class BlueprintBuilder;

class LayerHelper {
 public:
  using LayerType = Acts::Experimental::LayerBlueprintNode::LayerType;
  using Customizer =
      std::function<std::shared_ptr<Acts::Experimental::LayerBlueprintNode>(
          const dd4hep::DetElement&,
          std::shared_ptr<Acts::Experimental::LayerBlueprintNode>)>;

  explicit LayerHelper(const BlueprintBuilder& builder) : m_builder{&builder} {}

  LayerHelper& setLayerType(LayerType layerType) {
    m_layerType = layerType;
    return *this;
  }

  LayerHelper& endcap() { return setLayerType(LayerType::Disc); }
  LayerHelper& barrel() { return setLayerType(LayerType::Cylinder); }

  LayerHelper& setAxes(const std::string& axes) {
    m_axes = axes;
    return *this;
  }

  LayerHelper& setLayerAxes(const std::string& layerAxes) {
    m_layerAxes = layerAxes;
    return *this;
  }

  LayerHelper& setPattern(const std::string& pattern) {
    return setPattern(std::regex{pattern});
  }

  LayerHelper& setPattern(const std::regex& pattern) {
    m_pattern = pattern;
    return *this;
  }

  LayerHelper& setContainer(const dd4hep::DetElement& container) {
    m_container = container;
    return *this;
  }

  LayerHelper& setContainer(const std::string& name);

  LayerHelper& setEnvelope(const Acts::ExtentEnvelope& envelope) {
    m_envelope = envelope;
    return *this;
  }

  LayerHelper& setEmptyOk(bool emptyOk) {
    m_emptyOk = emptyOk;
    return *this;
  }

  LayerHelper& customize(Customizer customizer) {
    m_customizer = std::move(customizer);
    return *this;
  }

  auto& setAttachmentStrategy(
      std::optional<Acts::VolumeAttachmentStrategy> strategy) {
    m_attachmentStrategy = strategy;
    return *this;
  }

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> build()
      const;

  void addTo(Acts::Experimental::BlueprintNode& node) const;

 private:
  const ActsPlugins::DD4hep::BlueprintBuilder* m_builder;
  std::optional<LayerType> m_layerType;
  std::optional<std::string> m_axes;
  std::optional<std::string> m_layerAxes;
  std::optional<std::regex> m_pattern;
  std::optional<dd4hep::DetElement> m_container;
  std::optional<Acts::ExtentEnvelope> m_envelope;
  std::optional<Acts::VolumeAttachmentStrategy> m_attachmentStrategy;
  bool m_emptyOk = false;

  Customizer m_customizer;
};

class BarrelEndcapAssemblyHelper {
 public:
  explicit BarrelEndcapAssemblyHelper(const BlueprintBuilder& builder)
      : m_builder{&builder} {}

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> build()
      const;

  void addTo(Acts::Experimental::BlueprintNode& node) const;

  auto& customize(LayerHelper::Customizer customizer) {
    m_customizer = std::move(customizer);
    return *this;
  }

  auto& setAssembly(const dd4hep::DetElement& assembly) {
    m_assembly = assembly;
    return *this;
  }

  auto& setAxes(const std::string& barrel, const std::string& endcap) {
    m_barrelAxes = barrel;
    m_endcapAxes = endcap;
    return *this;
  }

  auto& setEndcapAxes(const std::string& axes) {
    m_endcapAxes = axes;
    return *this;
  }

  auto& setLayerPattern(const std::regex& pattern) {
    m_layerPattern = pattern;
    return *this;
  }

  auto& setAttachmentStrategies(Acts::VolumeAttachmentStrategy barrel,
                                Acts::VolumeAttachmentStrategy endcap) {
    m_barrelAttachmentStrategy = barrel;
    m_endcapAttachmentStrategy = endcap;
    return *this;
  }

 private:
  LayerHelper::Customizer m_customizer;

  std::optional<dd4hep::DetElement> m_assembly;
  std::optional<std::string> m_barrelAxes;
  std::optional<std::string> m_endcapAxes;
  std::optional<std::regex> m_layerPattern;
  std::optional<Acts::VolumeAttachmentStrategy> m_barrelAttachmentStrategy;
  std::optional<Acts::VolumeAttachmentStrategy> m_endcapAttachmentStrategy;
  const BlueprintBuilder* m_builder;
};

class BlueprintBuilder {
 public:
  using ElementFactory = std::function<std::shared_ptr<DD4hepDetectorElement>(
      const dd4hep::DetElement& detElement, const std::string& axes,
      double lengthScale)>;

  static std::shared_ptr<DD4hepDetectorElement> defaultElementFactory(
      const dd4hep::DetElement& detElement, const std::string& axes,
      double lengthScale);

  struct Config {
    ElementFactory elementFactory = defaultElementFactory;
    const dd4hep::Detector* dd4hepDetector;
    double lengthScale = 1.0;
  };

  explicit BlueprintBuilder(const Config& cfg,
                            std::unique_ptr<const Acts::Logger> logger_ =
                                Acts::getDefaultLogger("BlueprintBuilder",
                                                       Acts::Logging::INFO))
      : m_cfg(cfg), m_logger{std::move(logger_)} {
    if (m_cfg.dd4hepDetector == nullptr) {
      throw std::invalid_argument(
          "BlueprintBuilder: dd4hepDetector in config is null");
    }
  }

  std::shared_ptr<DD4hepDetectorElement> createDetectorElement(
      const dd4hep::DetElement& detElement, const std::string& axes) const;

  [[deprecated("Renamed to `makeLayer()`")]]
  std::shared_ptr<Acts::Experimental::LayerBlueprintNode> addLayer(
      const dd4hep::DetElement& detElement, const std::string& axes,
      std::optional<std::string> layerAxes = std::nullopt) {
    return makeLayer(detElement, axes, std::move(layerAxes));
  }

  std::shared_ptr<Acts::Experimental::LayerBlueprintNode> makeLayer(
      const dd4hep::DetElement& detElement, const std::string& axes,
      std::optional<std::string> layerAxes = std::nullopt) const;

  std::shared_ptr<Acts::Experimental::LayerBlueprintNode> makeLayer(
      const dd4hep::DetElement& parent,
      std::span<const dd4hep::DetElement> sensitives, const std::string& axes,
      std::optional<std::string> layerName = std::nullopt,
      std::optional<std::string> layerAxes = std::nullopt) const;

  std::shared_ptr<Acts::Experimental::StaticBlueprintNode> makeBeampipe() const;

  [[deprecated("Consider using .layerHelper() to produce the layers")]]
  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> addLayers(
      const dd4hep::DetElement& container, const std::string& axes,
      Acts::AxisDirection direction, const std::regex& layerPattern,
      const Acts::ExtentEnvelope& envelope = Acts::ExtentEnvelope::Zero());

  LayerHelper layerHelper() const;
  BarrelEndcapAssemblyHelper barrelEndcapAssemblyHelper() const;

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode>
  makeBarrelEndcapAssembly(
      const dd4hep::DetElement& assembly, const std::regex& layerPattern,
      const std::string& barrelAxes, const std::string& endcapAxes,
      const LayerHelper::Customizer& customizer = {}) const;

  static std::optional<dd4hep::DetElement> findDetElementByName(
      const dd4hep::DetElement& parent, const std::string& name);

  std::optional<dd4hep::DetElement> findDetElementByName(
      const std::string& name) const;

  std::string getPathToElementName(const dd4hep::DetElement& elem,
                                   const std::string& separator = "|") const;

  static std::vector<dd4hep::DetElement> findDetElementByNamePattern(
      const dd4hep::DetElement& parent, const std::regex& pattern);

  const Acts::Logger& logger() const { return *m_logger; }

 private:
  static std::vector<dd4hep::DetElement> resolveSensitives(
      const dd4hep::DetElement& detElement);

  dd4hep::DetElement world() const;

  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace DD4hep

}  // namespace ActsPlugins
