// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"

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
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceBoundsJsonConverter.hpp"

#include <memory>

namespace {

// @brief Get string representation of the surface bounds kind
//
// @return string representation of the surface bounds kind
template <typename bounds_t>
std::string getSurfaceBoundsKind() {
  if (std::is_same_v<bounds_t, Acts::EllipseBounds>) {
    return "Ellipse";
  } else if (std::is_same_v<bounds_t, Acts::RectangleBounds>) {
    return "Rectangle";
  } else if (std::is_same_v<bounds_t, Acts::TrapezoidBounds>) {
    return "Trapezoid";
  } else if (std::is_same_v<bounds_t, Acts::AnnulusBounds>) {
    return "Annulus";
  } else if (std::is_same_v<bounds_t, Acts::RadialBounds>) {
    return "Radial";
  } else if (std::is_same_v<bounds_t, Acts::DiscTrapezoidBounds>) {
    return "DiscTrapezoid";
  } else if (std::is_same_v<bounds_t, Acts::CylinderBounds>) {
    return "Cylinder";
  } else if (std::is_same_v<bounds_t, Acts::ConeBounds>) {
    return "Cone";
  } else if (std::is_same_v<bounds_t, Acts::LineBounds>) {
    return "Line";
  } else if (std::is_same_v<bounds_t, Acts::SurfaceBounds>) {
    return "DefaultBounds";
  } else if (std::is_same_v<bounds_t, Acts::InfiniteBounds>) {
    return "Infinite";
  } else {
    throw std::invalid_argument("Unknown surface bounds kind");
  }
}

// @brief Get string representation of the surface kind
//
// @return string representation of the surface kind
template <typename surface_t>
std::string getSurfaceKind() {
  if (std::is_same_v<surface_t, Acts::PlaneSurface>) {
    return "Plane";
  } else if (std::is_same_v<surface_t, Acts::DiscSurface>) {
    return "Disc";
  } else if (std::is_same_v<surface_t, Acts::CylinderSurface>) {
    return "Cylinder";
  } else if (std::is_same_v<surface_t, Acts::ConeSurface>) {
    return "Cone";
  } else if (std::is_same_v<surface_t, Acts::StrawSurface>) {
    return "Straw";
  } else if (std::is_same_v<surface_t, Acts::PerigeeSurface>) {
    return "Perigee";
  } else {
    throw std::invalid_argument("Unknown surface kind");
  }
}

// @brief Type-based surface bounds json encoding
//
// @tparam bounds_t surface bounds type
//
// @param bounds surface bounds to be converted
//
// @return json representation of the surface bounds
template <typename bounds_t>
nlohmann::json surfaceBoundsToJsonT(const bounds_t& bounds) {
  nlohmann::json jBounds = Acts::SurfaceBoundsJsonConverter::toJson(bounds);
  jBounds["kind"] = getSurfaceBoundsKind<bounds_t>();
  return jBounds;
}

// @brief Type-based surface json encoding
//
// @tparam surface_t surface type
//
// @param bounds surface to be converted
// @param gctx geometry context
// @param opt surface json conversion options
//
// @return json representation of the surface bounds
template <typename surface_t>
nlohmann::json surfaceToJsonT(const surface_t& surface,
                              const Acts::GeometryContext& gctx,
                              const Acts::SurfaceJsonConverter::Options& opt) {
  nlohmann::json jSurface;
  const auto sTransform = surface.localToGlobalTransform(gctx);

  jSurface["transform"] =
      Acts::Transform3JsonConverter::toJson(sTransform, opt.transformOptions);
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

// @brief Type-based surface json decoding
//
// @tparam surface_t surface type
//
// @param j json encoding of the surface
//
// @return shared pointer to the decoded surface
template <typename surface_t>
std::shared_ptr<Acts::Surface> surfaceFromJsonT(const nlohmann::json& j) {
  nlohmann::json jTransform = j["transform"];
  Acts::Transform3 sTransform =
      Acts::Transform3JsonConverter::fromJson(jTransform);
  return Acts::Surface::makeShared<surface_t>(sTransform);
}

// @brief Type-based surface json decoding
//
// @tparam surface_t surface type
// @tparam bounds_t surface bounds type
//
// @param j json encoding of the surface
//
// @return shared pointer to the decoded surface
template <typename surface_t, typename bounds_t>
std::shared_ptr<Acts::Surface> surfaceFromJsonT(const nlohmann::json& j) {
  nlohmann::json jTransform = j["transform"];
  Acts::Transform3 sTransform =
      Acts::Transform3JsonConverter::fromJson(jTransform);
  nlohmann::json jBounds = j["bounds"];
  auto sBounds = Acts::SurfaceBoundsJsonConverter::fromJson<bounds_t>(jBounds);
  return Acts::Surface::makeShared<surface_t>(sTransform, std::move(sBounds));
}

}  // namespace

void Acts::to_json(nlohmann::json& j,
                   const Acts::SurfaceAndMaterialWithContext& surface) {
  toJson(j, std::get<0>(surface), std::get<2>(surface));
  to_json(j, std::get<1>(surface).get());
}

void Acts::to_json(nlohmann::json& j, const Acts::Surface& surface) {
  Acts::GeometryContext gctx =
      Acts::GeometryContext::dangerouslyDefaultConstruct();
  j = SurfaceJsonConverter::toJson(gctx, surface);
}

void Acts::to_json(nlohmann::json& j,
                   const std::shared_ptr<const Acts::Surface>& surface) {
  Acts::GeometryContext gctx =
      Acts::GeometryContext::dangerouslyDefaultConstruct();
  j = SurfaceJsonConverter::toJson(gctx, *surface);
}

void Acts::toJson(nlohmann::json& j,
                  const std::shared_ptr<const Acts::Surface>& surface,
                  const Acts::GeometryContext& gctx) {
  j = SurfaceJsonConverter::toJson(gctx, *surface);
}

Acts::SurfaceJsonConverter::Config
Acts::SurfaceJsonConverter::Config::defaultConfig() {
  Config cfg;

  // Encoders
  cfg.surfaceEncoder.registerFunction(surfaceToJsonT<PlaneSurface>);
  cfg.surfaceEncoder.registerFunction(surfaceToJsonT<DiscSurface>);
  cfg.surfaceEncoder.registerFunction(surfaceToJsonT<CylinderSurface>);
  cfg.surfaceEncoder.registerFunction(surfaceToJsonT<ConeSurface>);
  cfg.surfaceEncoder.registerFunction(surfaceToJsonT<StrawSurface>);
  cfg.surfaceEncoder.registerFunction(surfaceToJsonT<PerigeeSurface>);

  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<EllipseBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<RectangleBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<TrapezoidBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<AnnulusBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(surfaceBoundsToJsonT<RadialBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<DiscTrapezoidBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<CylinderBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(surfaceBoundsToJsonT<ConeBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(surfaceBoundsToJsonT<LineBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<InfiniteBounds>);

  // Decoders
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<PlaneSurface>() + getSurfaceBoundsKind<EllipseBounds>(),
      surfaceFromJsonT<PlaneSurface, EllipseBounds>);
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<PlaneSurface>() + getSurfaceBoundsKind<RectangleBounds>(),
      surfaceFromJsonT<PlaneSurface, RectangleBounds>);
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<PlaneSurface>() + getSurfaceBoundsKind<TrapezoidBounds>(),
      surfaceFromJsonT<PlaneSurface, TrapezoidBounds>);
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<PlaneSurface>() + getSurfaceBoundsKind<InfiniteBounds>(),
      surfaceFromJsonT<PlaneSurface>);

  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<DiscSurface>() + getSurfaceBoundsKind<AnnulusBounds>(),
      surfaceFromJsonT<DiscSurface, AnnulusBounds>);
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<DiscSurface>() + getSurfaceBoundsKind<RadialBounds>(),
      surfaceFromJsonT<DiscSurface, RadialBounds>);
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<DiscSurface>() +
          getSurfaceBoundsKind<DiscTrapezoidBounds>(),
      surfaceFromJsonT<DiscSurface, DiscTrapezoidBounds>);

  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<CylinderSurface>() +
          getSurfaceBoundsKind<CylinderBounds>(),
      surfaceFromJsonT<CylinderSurface, CylinderBounds>);
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<ConeSurface>() + getSurfaceBoundsKind<ConeBounds>(),
      surfaceFromJsonT<ConeSurface, ConeBounds>);
  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<StrawSurface>() + getSurfaceBoundsKind<LineBounds>(),
      surfaceFromJsonT<StrawSurface, LineBounds>);

  cfg.surfaceDecoder.registerKind(
      getSurfaceKind<PerigeeSurface>() + getSurfaceBoundsKind<InfiniteBounds>(),
      surfaceFromJsonT<PerigeeSurface>);
  return cfg;
}

Acts::SurfaceJsonConverter::Config Acts::SurfaceJsonConverter::m_cfg =
    Acts::SurfaceJsonConverter::Config::defaultConfig();

std::shared_ptr<Acts::Surface> Acts::SurfaceJsonConverter::fromJson(
    const nlohmann::json& j) {
  std::shared_ptr<Acts::Surface> mutableSf = nullptr;
  mutableSf = m_cfg.surfaceDecoder(j);

  if (j.find("geo_id") != j.end() && !j["geo_id"].empty()) {
    GeometryIdentifier geoID = j["geo_id"].get<GeometryIdentifier>();
    mutableSf->assignGeometryId(geoID);
  } else {
    mutableSf->assignGeometryId(GeometryIdentifier(0));
  }
  mutableSf->assignIsSensitive(j["sensitive"].get<bool>());

  if (j.find("material") != j.end() && !j["material"].empty()) {
    const ISurfaceMaterial* surfaceMaterial = nullptr;
    from_json(j, surfaceMaterial);
    std::shared_ptr<const ISurfaceMaterial> sharedSurfaceMaterial(
        surfaceMaterial);
    mutableSf->assignSurfaceMaterial(sharedSurfaceMaterial);
  }
  return mutableSf;
}

nlohmann::json Acts::SurfaceJsonConverter::toJson(const GeometryContext& gctx,
                                                  const Surface& surface,
                                                  const Options& options) {
  nlohmann::json jSurface = m_cfg.surfaceEncoder(surface, gctx, options);
  nlohmann::json jBounds = m_cfg.surfaceBoundsEncoder(surface.bounds());
  jSurface["bounds"] = jBounds;
  jSurface["kind"] =
      jSurface["kind"].get<std::string>() + jBounds["kind"].get<std::string>();
  return jSurface;
}

nlohmann::json Acts::SurfaceJsonConverter::toJsonDetray(
    const GeometryContext& gctx, const Surface& surface,
    const Options& options) {
  nlohmann::json jSurface;
  const auto& sBounds = surface.bounds();
  const auto sTransform = surface.localToGlobalTransform(gctx);

  jSurface["transform"] =
      Transform3JsonConverter::toJson(sTransform, options.transformOptions);

  auto jMask =
      SurfaceBoundsJsonConverter::toJsonDetray(sBounds, options.portal);
  jSurface["mask"] = jMask;
  jSurface["source"] = surface.geometryId().value();
  jSurface["barcode"] = 0;
  jSurface["type"] =
      options.portal ? 0 : (surface.geometryId().sensitive() > 0 ? 1u : 2u);

  return jSurface;
}
