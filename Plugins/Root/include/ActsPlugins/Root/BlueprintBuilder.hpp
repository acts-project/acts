// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/TGeoAxes.hpp"
#include "ActsPlugins/Root/TGeoDetectorElement.hpp"

#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "TGeoMatrix.h"

class TGeoNode;

namespace ActsPlugins {

class TGeoBackend {
 public:
  static constexpr std::string_view kIdentifier = "TGeoBackend";

  struct NodeContext {
    const TGeoNode* node = nullptr;
    std::shared_ptr<const NodeContext> parent = nullptr;
  };

  struct Element {
    std::shared_ptr<const NodeContext> context = nullptr;

    friend bool operator==(const Element& lhs, const Element& rhs) {
      const NodeContext* lc = lhs.context.get();
      const NodeContext* rc = rhs.context.get();

      while (lc != nullptr && rc != nullptr) {
        if (lc->node != rc->node) {
          return false;
        }
        lc = lc->parent.get();
        rc = rc->parent.get();
      }

      return lc == nullptr && rc == nullptr;
    }
  };

  using AxisDefinition = TGeoAxes;
  struct LayerSpec {
    std::optional<AxisDefinition> axes;
    std::optional<AxisDefinition> layerAxes;
    std::optional<std::string> layerName;
  };

  using DetectorElement = TGeoDetectorElement;
  using DetectorElementPtr = std::shared_ptr<DetectorElement>;
  using DetectorElementFactory = std::function<DetectorElementPtr(
      const TGeoDetectorElement::Identifier&, const TGeoNode&,
      const TGeoMatrix&, AxisDefinition, double,
      std::shared_ptr<const Acts::ISurfaceMaterial>)>;
  using NameProvider = std::function<std::string(const Element&)>;
  using ElementPredicate = std::function<bool(const Element&)>;
  using IdentifierProvider =
      std::function<TGeoDetectorElement::Identifier(const Element&)>;

  static DetectorElementPtr defaultElementFactory(
      const TGeoDetectorElement::Identifier& identifier,
      const TGeoNode& tGeoNode, const TGeoMatrix& tGeoMatrix,
      AxisDefinition axes, double lengthScale,
      std::shared_ptr<const Acts::ISurfaceMaterial> material);

  struct Config {
    const TGeoNode* root = nullptr;
    double lengthScale = 1.0;
    DetectorElementFactory elementFactory = defaultElementFactory;
    NameProvider nameProvider = {};
    ElementPredicate sensitivePredicate = {};
    IdentifierProvider identifierProvider = {};
  };

  explicit TGeoBackend(const Config& cfg, const Acts::Logger& logger);

  DetectorElementPtr createDetectorElement(const Element& element,
                                           AxisDefinition axes) const;

  std::vector<std::shared_ptr<Acts::Surface>> makeSurfaces(
      std::span<const Element> sensitives, const LayerSpec& layerSpec) const;

  std::optional<Acts::Transform3> lookupLayerTransform(
      const Element& element, const LayerSpec& layerSpec) const;

  Element world() const;
  std::string nameOf(const Element& element) const;
  std::vector<Element> children(const Element& parent) const;
  Element parent(const Element& element) const;

  bool isSensitive(const Element& element) const;

  const TGeoNode& nodeOf(const Element& element) const;
  TGeoHMatrix transformOf(const Element& element) const;
  std::string pathOf(const Element& element) const;

  const Acts::Logger& logger() const { return *m_logger; }

 private:
  const NodeContext& contextOf(const Element& element) const;
  Element makeElement(const TGeoNode& node,
                      std::shared_ptr<const NodeContext> parent) const;
  TGeoDetectorElement::Identifier defaultIdentifierFor(
      const Element& element) const;

  Config m_cfg;
  const Acts::Logger* m_logger = nullptr;
  Element m_world = {};
  mutable std::vector<DetectorElementPtr> m_detectorElementStore = {};
};

using BlueprintBuilder = Acts::Experimental::BlueprintBuilder<TGeoBackend>;
using ElementLayerAssembler =
    Acts::Experimental::ElementLayerAssembler<TGeoBackend>;
using SensorLayerAssembler =
    Acts::Experimental::SensorLayerAssembler<TGeoBackend>;
using SensorLayer = Acts::Experimental::SensorLayer<TGeoBackend>;

}  // namespace ActsPlugins

namespace Acts::Experimental {
extern template class BlueprintBuilder<ActsPlugins::TGeoBackend>;
extern template class ElementLayerAssembler<ActsPlugins::TGeoBackend>;
extern template class SensorLayerAssembler<ActsPlugins::TGeoBackend>;
extern template class SensorLayer<ActsPlugins::TGeoBackend>;
}  // namespace Acts::Experimental
