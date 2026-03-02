// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/OpenDataDetector.hpp"

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
// Current (real) builder — will become BlueprintBuilder<DD4hepBackend> via alias
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
// Current TGeoAxes location — would move to ActsPlugins/GeometryAxes.hpp
#include "ActsPlugins/Root/TGeoAxes.hpp"
#include <Acts/Geometry/Blueprint.hpp>
#include <Acts/Geometry/BlueprintOptions.hpp>
#include <Acts/Geometry/BoundarySurfaceFace.hpp>
#include <Acts/Geometry/ContainerBlueprintNode.hpp>
#include <Acts/Geometry/CylinderVolumeBounds.hpp>
#include <Acts/Geometry/Extent.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Navigation/CylinderNavigationPolicy.hpp>
#include <Acts/Navigation/SurfaceArrayNavigationPolicy.hpp>
#include <Acts/Navigation/TryAllNavigationPolicy.hpp>
#include <Acts/Surfaces/SurfaceArray.hpp>
#include <Acts/Utilities/AxisDefinitions.hpp>
#include <ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp>

#include <algorithm>
#include <format>
#include <iterator>
#include <utility>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>

// ============================================================================
// MOCK — illustrating the proposed backend-abstraction design.
//
// In real code these declarations would live in:
//
//   ActsPlugins/GeometryAxes.hpp
//     — TGeoAxes renamed and lifted here (zero external deps, trivial move)
//
//   ActsPlugins/GeometryBackend.hpp
//     — shared concept; no geometry-source dependency
//
//   ActsPlugins/DD4hep/DD4hepBackend.hpp
//     — DD4hep-specific backend struct with LayerSpec
//
//   ActsPlugins/DD4hep/BlueprintBuilder.hpp  (updated)
//     — type aliases keep existing names unchanged for callers
// ============================================================================
namespace mock {

// Renamed from TGeoAxes. Same consteval construction, same validation.
// The name "TGeo" was always an accident of where the class first appeared —
// the concept (axis permutation encoding) is geometry-source-agnostic.
using GeometryAxes = ActsPlugins::TGeoAxes;

// Shared concept. Geometry-source backends must satisfy this to be usable
// with the generic LayerHelper<B> and BarrelEndcapAssemblyHelper<B>.
//
// Deliberately excludes makeLayer from the concept: the construction signature
// varies per backend (DD4hep/TGeo need axes, GeoModel derives them from shape),
// so it is enforced backend-internally rather than in the shared concept.
template <typename B>
concept GeometryBackend =
    requires {
      typename B::Element;   // native handle — value, pointer, or smart pointer
      typename B::LayerSpec; // backend-specific construction config
    } &&
    requires(const B& b, const typename B::Element& e, const std::string& name,
             const std::regex& pat) {
      // Tree traversal — name/pattern based, consistent across all backends
      { b.findByName(e, name) } -> std::same_as<std::optional<typename B::Element>>;
      { b.findByPattern(e, pat) } -> std::same_as<std::vector<typename B::Element>>;
      { b.nameOf(e) } -> std::convertible_to<std::string_view>;
      // Role detection — DD4hep uses type flags; TGeo/GeoModel use naming
      { b.isBarrel(e) } -> std::convertible_to<bool>;
      { b.isEndcap(e) } -> std::convertible_to<bool>;
      // Root of the world tree, for convenience findByName(name) overloads
      { b.world() } -> std::same_as<typename B::Element>;
    };

// DD4hep backend. Moves the construction knowledge that currently lives
// inside BlueprintBuilder out into an explicit, swappable struct.
struct DD4hepBackend {
  using Element = dd4hep::DetElement;

  // LayerSpec carries the axes configuration needed to construct a detector
  // element from a DD4hep node. TGeoBackend::LayerSpec would be identical
  // (both use TGeo internally). GeoModelBackend::LayerSpec would be empty —
  // its converters derive orientation from the shape geometry automatically.
  struct LayerSpec {
    GeometryAxes axes;
    std::optional<GeometryAxes> layerAxes;
  };

  // Config: moves here from BlueprintBuilder::Config (same fields)
  // Methods: findByName, findByPattern, nameOf, isBarrel, isEndcap, world,
  //          makeLayer(Element, LayerSpec), makeBeampipe()
  // isBarrel/isEndcap implemented using DD4hep type flags (reliable, as
  // discussed — not name-based)
};

// Type aliases live in ActsPlugins/DD4hep/BlueprintBuilder.hpp.
// Existing callers keep the unqualified names unchanged — zero breakage.
//
//   using DD4hepBlueprintBuilder =
//       ActsPlugins::BlueprintBuilder<DD4hepBackend>;
//   using DD4hepLayerHelper =
//       ActsPlugins::LayerHelper<DD4hepBackend>;
//   using DD4hepBarrelEndcapAssemblyHelper =
//       ActsPlugins::BarrelEndcapAssemblyHelper<DD4hepBackend>;
//
// A TGeo equivalent would look like:
//
//   struct TGeoBackend {
//     using Element   = TGeoNode*;
//     using LayerSpec = DD4hepBackend::LayerSpec; // identical — same axes concept
//     // isBarrel/isEndcap: name-pattern based (configurable)
//   };
//   using TGeoBlueprintBuilder = ActsPlugins::BlueprintBuilder<TGeoBackend>;
//
// And GeoModel:
//
//   struct GeoModelBackend {
//     using Element   = PVConstLink;
//     struct LayerSpec {};  // empty — orientation derived from shape
//     // isBarrel/isEndcap: TBD, likely name-based initially
//   };
//   using GeoModelBlueprintBuilder = ActsPlugins::BlueprintBuilder<GeoModelBackend>;

}  // namespace mock
// ============================================================================
// END MOCK
// ============================================================================

namespace ActsExamples {
OpenDataDetector::OpenDataDetector(const Config& cfg,
                                   const Acts::GeometryContext& gctx)
    : DD4hepDetectorBase{cfg}, m_cfg{cfg} {
  ACTS_INFO("OpenDataDetector construct");
  construct(gctx);
}

auto OpenDataDetector::config() const -> const Config& {
  return m_cfg;
}

std::shared_ptr<ActsPlugins::DD4hepDetectorElement>
OpenDataDetector::defaultDetectorElementFactory(
    const dd4hep::DetElement& element, ActsPlugins::TGeoAxes axes,
    double scale) {
  return std::make_shared<ActsPlugins::DD4hepDetectorElement>(element, axes,
                                                              scale);
}

namespace {

int parseLayerNumber(const dd4hep::DetElement& elem,
                     const std::regex& pattern) {
  std::cmatch match;

  if (!std::regex_match(elem.name(), match, pattern)) {
    throw std::runtime_error(std::format(
        "Layer name {} does not match expected pattern", elem.name()));
  }

  if (match.size() != 2) {
    throw std::runtime_error(std::format(
        "Layer name {} matched pattern but did not capture layer number",
        elem.name()));
  }

  int n = std::stoi(match[1]);
  return n;
}

}  // namespace

void OpenDataDetector::construct(const Acts::GeometryContext& gctx) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using namespace Acts::UnitLiterals;
  using enum AxisDirection;

  // In the new design this would be DD4hepBlueprintBuilder — a type alias for
  // BlueprintBuilder<DD4hepBackend>. Config structure is unchanged.
  ActsPlugins::DD4hep::BlueprintBuilder builder{
      {
          .dd4hepDetector = &dd4hepDetector(),
          .lengthScale = Acts::UnitConstants::cm,
          .gctx = gctx,
      },
      logger().cloneWithSuffix("BlpBld")};

  // BARREL: XYZ
  // ENDCAP: XZY

  Blueprint::Config cfg;
  cfg.envelope[AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisR] = {0_mm, 20_mm};
  Blueprint root{cfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);

  outer.addChild(builder.makeBeampipe());

  using AttachmentStrategy = Acts::VolumeAttachmentStrategy;
  using ResizeStrategy = Acts::VolumeResizeStrategy;
  using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;

  auto constant = [this]<typename... Args>(std::format_string<Args...> fmt,
                                           Args&&... values) -> int {
    return dd4hepDetector().constant<int>(
        std::format(fmt, std::forward<Args>(values)...));
  };

  // Customizer type would be DD4hepLayerHelper::Customizer (via alias).
  // The lambda body is unchanged — it receives the native dd4hep::DetElement,
  // so experiment code retains full access to backend-specific data.
  // Generic lambdas (auto& elem) would also work for cross-backend customizers
  // that only operate on ACTS node types.
  auto makeCustomizer = [&](const std::string& det,
                            const std::regex& layerPattern)
      -> ActsPlugins::DD4hep::LayerHelper::Customizer {
    return [&, det](const auto& elem, auto layer) {
      int n = parseLayerNumber(elem, layerPattern);

      std::string name = elem.name();
      using enum SrfArrayNavPol::LayerType;
      SrfArrayNavPol::Config cfg;

      bool isEndcap = name.find("Endcap") != std::string::npos;

      if (isEndcap) {
        cfg.bins = {constant("{}_e_sf_b_r", det),
                    constant("{}_e_sf_b_phi", det)};
        cfg.layerType = Disc;
      } else {
        cfg.bins = {constant("{}_b{}_sf_b_phi", det, n),
                    constant("{}_b_sf_b_z", det)};
        cfg.layerType = Cylinder;
      }

      layer->setNavigationPolicyFactory(NavigationPolicyFactory{}
                                            .add<CylinderNavigationPolicy>()
                                            .add<SrfArrayNavPol>(cfg)
                                            .asUniquePtr());

      layer->setEnvelope(
          Acts::ExtentEnvelope::Zero().set(AxisZ, {2, 2}).set(AxisR, {2, 2}));

      return layer;
    };
  };

  // Generic lambda — auto& elem means this works with any backend's
  // BarrelEndcapAssemblyHelper::Customizer without modification.
  auto customizeContainer = [&](const auto& /*elem*/, auto node) {
    node->setAttachmentStrategy(AttachmentStrategy::Gap);
    node->setResizeStrategies(ResizeStrategy::Gap, ResizeStrategy::Gap);
    return node;
  };

  // Key API changes vs. today:
  //
  //   setAssembly("Pixels")         — string overload added; no need to call
  //                                   findDetElementByName externally first
  //
  //   setLayerFilter(layerPattern)  — renamed from setLayerPattern; name now
  //                                   reflects that only matching is used here,
  //                                   not capture groups (those remain the
  //                                   caller's responsibility in customizeLayer)
  //
  //   setAxes("XYZ", "XZY")        — unchanged, but now conditionally compiled
  //                                   via `requires` on Backend::LayerSpec
  //                                   having an axes field; GeoModel callers
  //                                   simply don't call this method

  std::regex pixelLayerPattern{"(?:PixelLayer|PixelEndcap[NP])(\\d)"};
  builder.barrelEndcapAssemblyHelper()
      .setAssembly("Pixels")
      .setAxes("XYZ", "XZY")
      .setLayerFilter(pixelLayerPattern)
      .customizeLayer(makeCustomizer("pix", pixelLayerPattern))
      .customize(customizeContainer)
      .addTo(outer);

  std::regex sstripLayerPattern{
      "(?:ShortStripLayer|ShortStripEndcap[NP])(\\d)"};
  builder.barrelEndcapAssemblyHelper()
      .setAssembly("ShortStrips")
      .setAxes("XYZ", "XZY")
      .setLayerFilter(sstripLayerPattern)
      .customizeLayer(makeCustomizer("ss", sstripLayerPattern))
      .customize(customizeContainer)
      .addTo(outer);

  std::regex lstripLayerPattern{"(?:LongStripLayer|LongStripEndcap[NP])(\\d)"};
  builder.barrelEndcapAssemblyHelper()
      .setAssembly("LongStrips")
      .setAxes("XYZ", "XZY")
      .setLayerFilter(lstripLayerPattern)
      .customizeLayer(makeCustomizer("ls", lstripLayerPattern))
      .customize(customizeContainer)
      .addTo(outer);

  BlueprintOptions options;

  m_trackingGeometry = root.construct(options, gctx, logger());
}

}  // namespace ActsExamples
