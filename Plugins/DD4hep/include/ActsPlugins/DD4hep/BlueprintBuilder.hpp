// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/DD4hep/LayerBlueprintNode.hpp"

#include <functional>
#include <memory>
#include <regex>
#include <string>

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

  std::shared_ptr<LayerBlueprintNode> addLayer(
      const dd4hep::DetElement& detElement, const std::string& axes);

  std::shared_ptr<Acts::Experimental::CylinderContainerBlueprintNode> addLayers(
      const dd4hep::DetElement& container, const std::string& axes,
      Acts::AxisDirection direction, const std::regex& layerPattern);

  static std::optional<dd4hep::DetElement> findDetElementByName(
      const dd4hep::DetElement& parent, const std::string& name);

  std::optional<dd4hep::DetElement> findDetElementByName(
      const std::string& name);

  static std::vector<dd4hep::DetElement> findDetElementByNamePattern(
      const dd4hep::DetElement& parent, const std::regex& pattern);

 private:
  static std::vector<dd4hep::DetElement> resolveSensitives(
      const dd4hep::DetElement& detElement);

  dd4hep::DetElement world() const;

  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace DD4hep

}  // namespace ActsPlugins
