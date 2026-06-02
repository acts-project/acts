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
#include <string>
#include <string_view>
#include <vector>

#include "TGeoMatrix.h"

class TGeoNode;

namespace ActsPlugins {

/// Backend adapter for building ACTS blueprints from a TGeo geometry tree.
class TGeoBlueprintBuilderBackend {
 public:
  /// Identifier used by the generic blueprint infrastructure.
  static constexpr std::string_view kIdentifier = "TGeoBlueprintBuilderBackend";

  /// Context for one node in the TGeo hierarchy.
  struct NodeContext {
    /// Current TGeo node.
    const TGeoNode* node = nullptr;
    /// Parent node context in the hierarchy (null for the root element).
    std::shared_ptr<const NodeContext> parent = nullptr;
  };

  /// Lightweight backend element that points to a node context.
  struct Element {
    /// Shared node context for this element.
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

  /// Axis definitions used for detector element and layer surface creation.
  using AxisDefinition = TGeoAxes;
  /// Optional layer configuration used during layer assembly.
  struct LayerSpec {
    /// Optional override for sensitive element axes.
    std::optional<AxisDefinition> axes;
    /// Optional override for the assembled layer surface axes.
    std::optional<AxisDefinition> layerAxes;
    /// Optional override for the assembled layer name.
    std::optional<std::string> layerName;
  };

  /// Concrete detector element type used by this backend.
  using DetectorElement = TGeoDetectorElement;
  /// Shared pointer to the detector element type.
  using DetectorElementPtr = std::shared_ptr<DetectorElement>;
  /// Factory function used to instantiate detector elements.
  using DetectorElementFactory = std::function<DetectorElementPtr(
      const TGeoDetectorElement::Identifier&, const TGeoNode&,
      const TGeoMatrix&, AxisDefinition, double,
      std::shared_ptr<const Acts::ISurfaceMaterial>)>;
  /// Predicate used to classify geometry elements as sensitive.
  using ElementPredicate = std::function<bool(const Element&)>;
  /// Provider that maps an element to a detector element identifier.
  using IdentifierProvider =
      std::function<TGeoDetectorElement::Identifier(const Element&)>;

  /// Default detector element factory implementation.
  /// @param identifier Identifier assigned to the detector element.
  /// @param tGeoNode Source TGeo node.
  /// @param tGeoMatrix Global transform matrix of the source node.
  /// @param axes Axis definition used to build the ACTS surface.
  /// @param lengthScale Scale factor converting TGeo lengths to ACTS units.
  /// @param material Optional surface material assigned to created surfaces.
  /// @return Shared detector element instance.
  static DetectorElementPtr defaultElementFactory(
      const TGeoDetectorElement::Identifier& identifier,
      const TGeoNode& tGeoNode, const TGeoMatrix& tGeoMatrix,
      AxisDefinition axes, double lengthScale,
      std::shared_ptr<const Acts::ISurfaceMaterial> material);

  /// Runtime configuration for the TGeo blueprint backend.
  struct Config {
    /// Root node of the input TGeo geometry tree.
    const TGeoNode* root = nullptr;
    /// Scale factor converting TGeo lengths into ACTS units.
    double lengthScale = 1.0;
    /// Factory used to construct detector elements for sensitive nodes.
    DetectorElementFactory elementFactory = defaultElementFactory;
    /// Predicate deciding whether an element should be treated as sensitive.
    ElementPredicate sensitivePredicate = {};
    /// Provider used to compute detector element identifiers.
    IdentifierProvider identifierProvider = {};
  };

  /// Construct the backend from configuration and logger.
  /// @param cfg Backend configuration.
  /// @param logger Logger instance used for diagnostics.
  explicit TGeoBlueprintBuilderBackend(const Config& cfg,
                                       const Acts::Logger& logger);

  /// Create a detector element for a sensitive geometry element.
  /// @param element Sensitive element to convert.
  /// @param axes Axis definition used to build the detector element surface.
  /// @return Shared detector element instance.
  DetectorElementPtr createDetectorElement(const Element& element,
                                           AxisDefinition axes) const;

  /// Build ACTS surfaces for a set of sensitive elements.
  /// @param sensitives Sensitive elements that contribute surfaces.
  /// @param layerSpec Optional layer assembly overrides.
  /// @return Created ACTS surfaces.
  std::vector<std::shared_ptr<Acts::Surface>> makeSurfaces(
      std::span<const Element> sensitives, const LayerSpec& layerSpec) const;

  /// Look up an optional transform to use for the assembled layer.
  /// @param element Reference element for the lookup.
  /// @param layerSpec Optional layer assembly overrides.
  /// @return Layer transform if one can be determined.
  std::optional<Acts::Transform3> lookupLayerTransform(
      const Element& element, const LayerSpec& layerSpec) const;

  /// Access the world/root element of the configured TGeo tree.
  /// @return Root element of the geometry hierarchy.
  Element world() const;
  /// Retrieve the node name for an element.
  /// @param element Element whose node name is requested.
  /// @return Node name.
  std::string nameOf(const Element& element) const;
  /// Retrieve the direct children of a parent element.
  /// @param parent Parent element.
  /// @return Direct child elements.
  std::vector<Element> children(const Element& parent) const;
  /// Retrieve the parent element of an element.
  /// @param element Element whose parent is requested.
  /// @return Parent element.
  Element parent(const Element& element) const;

  /// Check whether an element is considered sensitive.
  /// @param element Element to classify.
  /// @return True if the element is sensitive.
  bool isSensitive(const Element& element) const;

  /// Access the underlying TGeo node for an element.
  /// @param element Element whose node is requested.
  /// @return Referenced TGeo node.
  const TGeoNode& nodeOf(const Element& element) const;
  /// Compute the world transform of an element.
  /// @param element Element whose transform is requested.
  /// @return Global transformation matrix.
  TGeoHMatrix transformOf(const Element& element) const;
  /// Compute the full geometry path of an element.
  /// @param element Element whose path is requested.
  /// @return Full geometry path.
  std::string pathOf(const Element& element) const;

  /// Access the backend logger.
  /// @return Backend logger instance.
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

/// Blueprint builder specialization for the TGeo backend.
using BlueprintBuilder =
    Acts::Experimental::BlueprintBuilder<TGeoBlueprintBuilderBackend>;
/// Element-layer assembler specialization for the TGeo backend.
using ElementLayerAssembler =
    Acts::Experimental::ElementLayerAssembler<TGeoBlueprintBuilderBackend>;
/// Sensor-layer assembler specialization for the TGeo backend.
using SensorLayerAssembler =
    Acts::Experimental::SensorLayerAssembler<TGeoBlueprintBuilderBackend>;
/// Sensor-layer alias for the TGeo backend.
using SensorLayer =
    Acts::Experimental::SensorLayer<TGeoBlueprintBuilderBackend>;

}  // namespace ActsPlugins

namespace Acts::Experimental {
extern template class BlueprintBuilder<
    ActsPlugins::TGeoBlueprintBuilderBackend>;
extern template class ElementLayerAssembler<
    ActsPlugins::TGeoBlueprintBuilderBackend>;
extern template class SensorLayerAssembler<
    ActsPlugins::TGeoBlueprintBuilderBackend>;
extern template class SensorLayer<ActsPlugins::TGeoBlueprintBuilderBackend>;
}  // namespace Acts::Experimental
