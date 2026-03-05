// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/TypeDispatcher.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/AlgebraJsonConverter.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"
#include "ActsPlugins/Json/JsonKindDispatcher.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceBoundsJsonConverter.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <Eigen/src/Core/util/Meta.h>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{
class ISurfaceMaterial;

using SurfaceAndMaterialWithContext =
    std::tuple<std::shared_ptr<const Surface>,
               std::shared_ptr<const ISurfaceMaterial>, GeometryContext>;

/// Convert SurfaceAndMaterialWithContext to JSON
/// @param j Destination JSON object
/// @param surface Source SurfaceAndMaterialWithContext to convert
void to_json(nlohmann::json& j, const SurfaceAndMaterialWithContext& surface);

/// Convert Surface to JSON
/// @param j Destination JSON object
/// @param surface Source Surface to convert
/// @note it will take the default context
void to_json(nlohmann::json& j, const Surface& surface);

/// Convert shared_ptr<Surface> to JSON
/// @param j Destination JSON object
/// @param surface Source shared_ptr<Surface> to convert
/// @note it will take the default context
void to_json(nlohmann::json& j, const std::shared_ptr<const Surface>& surface);

/// Contextual conversion of a surface
///
/// @param j the json to be filled
/// @param surface the surface to be converted
/// @param gctx the geometry context for this
void toJson(nlohmann::json& j, const std::shared_ptr<const Surface>& surface,
            const Acts::GeometryContext& gctx);

class SurfaceJsonConverter {
 public:
  using SurfaceBoundsEncoder = TypeDispatcher<SurfaceBounds, nlohmann::json()>;
  using SurfacePlanarBoundsEncoder =
      TypeDispatcher<PlanarBounds, nlohmann::json()>;
  using SurfaceDiscBoundsEncoder = TypeDispatcher<DiscBounds, nlohmann::json()>;

  struct Options {
    /// Transform serialization options
    Transform3JsonConverter::Options transformOptions =
        Transform3JsonConverter::Options{};
    /// Write material information
    bool writeMaterial = true;
    /// Write surface as portal
    bool portal = false;
  };

  using SurfaceEncoder =
      TypeDispatcher<Surface,
                     nlohmann::json(const GeometryContext&, const Options&)>;
  using SurfaceDecoder = JsonKindDispatcher<std::shared_ptr<Surface>>;

  struct Config {
    SurfaceEncoder surfaceEncoder{};
    SurfaceBoundsEncoder surfaceBoundsEncoder{};

    SurfaceDecoder surfaceDecoder{};

    static Config defaultConfig();
  };

  /// Contextual conversion of a surface
  ///
  /// @param gctx the geometry context for this
  /// @param surface the surface to be converted
  /// @param options the writing options for the surfaces
  ///
  /// @return a json object representing the surface
  static nlohmann::json toJson(
      const GeometryContext& gctx, const Surface& surface,
      const Options& options = Options{
          .transformOptions = Transform3JsonConverter::Options{},
          .writeMaterial = true,
          .portal = false});

  /// Contextual conversion of a surface - Detray export
  ///
  /// @param gctx the geometry context for this
  /// @param surface the surface to be converted
  /// @param options the writing options for the surfaces
  ///
  /// @note reading back detray json is not supported and will fail
  ///
  /// @return a json object representing the surface
  static nlohmann::json toJsonDetray(
      const GeometryContext& gctx, const Surface& surface,
      const Options& options = Options{
          .transformOptions = Transform3JsonConverter::Options{},
          .writeMaterial = true,
          .portal = false});

  /// @brief The Surface converter from json
  ///
  /// @param jSurface the surface json object
  ///
  /// @return a shared object created from json input
  static std::shared_ptr<Surface> fromJson(const nlohmann::json& jSurface);

  static Config m_cfg;
};

// This macro create a conversion for the surface type
NLOHMANN_JSON_SERIALIZE_ENUM(
    Surface::SurfaceType,
    {{Surface::SurfaceType::Cone, "ConeSurface"},
     {Surface::SurfaceType::Cylinder, "CylinderSurface"},
     {Surface::SurfaceType::Disc, "DiscSurface"},
     {Surface::SurfaceType::Perigee, "PerigeeSurface"},
     {Surface::SurfaceType::Plane, "PlaneSurface"},
     {Surface::SurfaceType::Straw, "StrawSurface"},
     {Surface::SurfaceType::Curvilinear, "CurvilinearSurface"},
     {Surface::SurfaceType::Other, "Other"}})

template <typename bounds_t>
std::string getSurfaceBoundsKind() {
  if (std::is_same_v<bounds_t, EllipseBounds>) {
    return "Ellipse";
  } else if (std::is_same_v<bounds_t, RectangleBounds>) {
    return "Rectangle";
  } else if (std::is_same_v<bounds_t, TrapezoidBounds>) {
    return "Trapezoid";
  } else if (std::is_same_v<bounds_t, AnnulusBounds>) {
    return "Annulus";
  } else if (std::is_same_v<bounds_t, RadialBounds>) {
    return "Radial";
  } else if (std::is_same_v<bounds_t, DiscTrapezoidBounds>) {
    return "DiscTrapezoid";
  } else if (std::is_same_v<bounds_t, CylinderBounds>) {
    return "Cylinder";
  } else if (std::is_same_v<bounds_t, ConeBounds>) {
    return "Cone";
  } else if (std::is_same_v<bounds_t, LineBounds>) {
    return "Line";
  } else if (std::is_same_v<bounds_t, SurfaceBounds>) {
    return "DefaultBounds";
  } else if (std::is_same_v<bounds_t, InfiniteBounds>) {
    return "Infinite";
  } else {
    throw std::invalid_argument("Unknown surface bounds kind");
  }
}

template <typename surface_t>
std::string getSurfaceKind() {
  if (std::is_same_v<surface_t, PlaneSurface>) {
    return "Plane";
  } else if (std::is_same_v<surface_t, DiscSurface>) {
    return "Disc";
  } else if (std::is_same_v<surface_t, CylinderSurface>) {
    return "Cylinder";
  } else if (std::is_same_v<surface_t, ConeSurface>) {
    return "Cone";
  } else if (std::is_same_v<surface_t, StrawSurface>) {
    return "Straw";
  } else if (std::is_same_v<surface_t, PerigeeSurface>) {
    return "Perigee";
  } else {
    throw std::invalid_argument("Unknown surface kind");
  }
}

template <typename bounds_t>
nlohmann::json surfaceBoundsToJsonT(const bounds_t& bounds) {
  nlohmann::json jBounds = SurfaceBoundsJsonConverter::toJson(bounds);
  jBounds["kind"] = getSurfaceBoundsKind<bounds_t>();
  return jBounds;
}

template <typename surface_t>
nlohmann::json surfaceToJsonT(const surface_t& surface,
                              const GeometryContext& gctx,
                              const SurfaceJsonConverter::Options& opt) {
  nlohmann::json jSurface;
  const auto sTransform = surface.localToGlobalTransform(gctx);

  jSurface["transform"] =
      Transform3JsonConverter::toJson(sTransform, opt.transformOptions);
  jSurface["type"] = surface.type();
  jSurface["geo_id"] = nlohmann::json(surface.geometryId());
  jSurface["sensitive"] = surface.isSensitive();
  if (surface.surfaceMaterial() != nullptr && opt.writeMaterial) {
    jSurface["material"] =
        nlohmann::json(surface.surfaceMaterial())["material"];
  }
  jSurface["kind"] = getSurfaceKind<surface_t>();
  return jSurface;
}

template <typename surface_t, typename bounds_t>
std::shared_ptr<Surface> surfaceFromJsonT(const nlohmann::json& j) {
  nlohmann::json jTransform = j["transform"];
  Acts::Transform3 sTransform =
      Acts::Transform3JsonConverter::fromJson(jTransform);
  nlohmann::json jBounds = j["bounds"];
  auto sBounds = SurfaceBoundsJsonConverter::fromJson<bounds_t>(jBounds);
  return Surface::makeShared<surface_t>(sTransform, std::move(sBounds));
}

/// @}
}  // namespace Acts
