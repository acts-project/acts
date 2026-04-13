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
#include <string_view>
#include <vector>

#include <DD4hep/Detector.h>

namespace dd4hep {
class DetElement;
}  // namespace dd4hep

namespace ActsPlugins::DD4hep {

/// Backend adapter that maps a DD4hep detector hierarchy onto the ACTS
/// experimental blueprint-builder interface.
///
/// The backend exposes DD4hep detector elements as blueprint input elements
/// and provides the conversions needed to create ACTS surfaces, detector
/// elements, and optional layer transforms.
class DD4hepBackend {
 public:
  /// Identifier string used in diagnostics emitted by the generic blueprint
  /// builder.
  static constexpr std::string_view kIdentifier = "DD4hepBackend";

  /// DD4hep detector-element handle type consumed by the generic builder.
  using Element = dd4hep::DetElement;
  /// Axis-definition type forwarded to ROOT-based surface conversion helpers.
  using AxisDefinition = TGeoAxes;
  /// Layer-specific configuration forwarded from the generic builder to the
  /// DD4hep backend.
  struct LayerSpec {
    /// Sensitive-surface axis convention.
    std::optional<AxisDefinition> axes;
    /// Optional axis convention used to derive the layer transform.
    std::optional<AxisDefinition> layerAxes;
    /// Optional explicit layer name override.
    std::optional<std::string> layerName;
  };
  /// Concrete ACTS detector-element implementation used for DD4hep geometry.
  using DetectorElement = DD4hepDetectorElement;
  /// Shared pointer to a DD4hep-backed ACTS detector element.
  using DetectorElementPtr = std::shared_ptr<DetectorElement>;
  /// Factory that creates detector elements for converted DD4hep sensitives.
  using DetectorElementFactory = std::function<DetectorElementPtr(
      const Element& detElement, AxisDefinition axes, double lengthScale)>;

  /// Default detector-element factory used when @ref Config::elementFactory is
  /// not overridden.
  /// @param detElement DD4hep sensitive detector element to wrap.
  /// @param axes Axis convention used for surface conversion.
  /// @param lengthScale Unit scale applied during geometry conversion.
  /// @return Newly created DD4hep-backed ACTS detector element.
  static DetectorElementPtr defaultElementFactory(const Element& detElement,
                                                  AxisDefinition axes,
                                                  double lengthScale);

  /// Configuration of the DD4hep backend instance.
  struct Config {
    /// Factory used to create ACTS detector elements from DD4hep sensitives.
    DetectorElementFactory elementFactory = defaultElementFactory;
    /// DD4hep detector description that owns the world hierarchy.
    const dd4hep::Detector* dd4hepDetector;
    /// Unit scale applied when converting DD4hep lengths to ACTS units.
    double lengthScale = 1.0;
    /// Geometry context used when constructing ACTS detector elements.
    std::reference_wrapper<const Acts::GeometryContext> gctx;
  };

  /// Construct the DD4hep backend.
  /// @param cfg Backend configuration and DD4hep detector handle.
  /// @param logger Logger used for diagnostics.
  explicit DD4hepBackend(const Config& cfg, const Acts::Logger& logger);

  /// Create an ACTS detector element from a DD4hep sensitive element.
  /// @param detElement DD4hep sensitive element to convert.
  /// @param axes Axis convention used for the converted surface.
  /// @return Shared pointer to the created detector element.
  DetectorElementPtr createDetectorElement(const Element& detElement,
                                           AxisDefinition axes) const;

  /// Convert a set of DD4hep sensitive elements into ACTS surfaces.
  /// @param sensitives Sensitive DD4hep elements belonging to one layer.
  /// @param layerSpec Layer configuration controlling axes and naming.
  /// @return Converted ACTS surfaces for the given sensitives.
  std::vector<std::shared_ptr<Acts::Surface>> makeSurfaces(
      std::span<const Element> sensitives, const LayerSpec& layerSpec) const;

  /// Derive the layer transform from a DD4hep detector element when possible.
  /// @param element DD4hep element providing the geometric context.
  /// @param layerSpec Layer configuration controlling the transform lookup.
  /// @return The derived layer transform, or `std::nullopt` if none is
  ///         available.
  std::optional<Acts::Transform3> lookupLayerTransform(
      const Element& element, const LayerSpec& layerSpec) const;

  /// Create a static beampipe blueprint node from the DD4hep world geometry.
  /// @return Shared pointer to the generated beampipe node.
  std::shared_ptr<Acts::Experimental::StaticBlueprintNode> makeBeampipe() const;

  /// Return the DD4hep world detector element.
  /// @return Root detector element of the DD4hep hierarchy.
  Element world() const;
  /// Return the fully qualified DD4hep name of an element.
  /// @param element DD4hep detector element to inspect.
  /// @return Full element name.
  std::string nameOf(const Element& element) const;
  /// Return the direct DD4hep child elements of a parent element.
  /// @param parent Parent DD4hep detector element.
  /// @return Direct children of @p parent.
  std::vector<Element> children(const Element& parent) const;
  /// Return the direct parent DD4hep element of a child element.
  /// @param element Child DD4hep detector element.
  /// @return Parent detector element.
  Element parent(const Element& element) const;

  /// Check whether a DD4hep element represents a sensitive detector element.
  /// @param element DD4hep detector element to classify.
  /// @return `true` if the element is sensitive.
  bool isSensitive(const Element& element) const;
  /// Check whether a DD4hep element represents a barrel sub-detector.
  /// @param element DD4hep detector element to classify.
  /// @return `true` if the element is tagged as barrel.
  bool isBarrel(const Element& element) const;
  /// Check whether a DD4hep element represents an endcap sub-detector.
  /// @param element DD4hep detector element to classify.
  /// @return `true` if the element is tagged as endcap.
  bool isEndcap(const Element& element) const;
  /// Check whether a DD4hep element is part of the tracking detector.
  /// @param element DD4hep detector element to classify.
  /// @return `true` if the element belongs to the tracker.
  bool isTracker(const Element& element) const;

  /// Retrieves a named integer constant from the DD4hep detector description.
  /// The name is constructed by formatting @p fmt with @p args.
  /// @tparam Args Types used to format the constant name.
  /// @param fmt Format string used to construct the DD4hep constant name.
  /// @param args Format arguments substituted into @p fmt.
  /// @return Integer constant value stored in the DD4hep detector description.
  template <typename... Args>
  int constant(std::format_string<Args...> fmt, Args&&... args) const {
    return m_cfg.dd4hepDetector->constant<int>(
        std::format(fmt, std::forward<Args>(args)...));
  }

  /// Return the logger associated with this backend.
  /// @return Logger used for diagnostics.
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
