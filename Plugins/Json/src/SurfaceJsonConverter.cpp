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
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/PointBounds.hpp"
#include "Acts/Surfaces/PointSurface.hpp"
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
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// @brief Registration key for a surface, from the serialized SurfaceType
std::string surfaceKey(Acts::Surface::SurfaceType type) {
  return nlohmann::json(type).get<std::string>();
}

// @brief Registration key for surface bounds, from the serialized BoundsType
std::string boundsKey(Acts::SurfaceBounds::BoundsType type) {
  return nlohmann::json(type).get<std::string>();
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
  return Acts::SurfaceBoundsJsonConverter::toJson(bounds);
}

// @brief Type-based surface json encoding
//
// @tparam surface_t surface type
//
// @param surface surface to be converted
// @param gctx geometry context
// @param opt surface json conversion options
//
// @return json representation of the surface
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
  if (surface.hasMaterial() && opt.writeMaterial) {
    jSurface["material"] =
        nlohmann::json(surface.surfaceMaterial())["material"];
  }
  return jSurface;
}

// @brief Type-based surface bounds json decoding
//
// @tparam bounds_t surface bounds type
//
// @param j json encoding of the surface bounds
//
// @return shared pointer to the decoded bounds, as the base type
template <typename bounds_t>
std::shared_ptr<const Acts::SurfaceBounds> boundsFromJsonT(
    const nlohmann::json& j) {
  return Acts::SurfaceBoundsJsonConverter::fromJson<bounds_t>(j);
}

// @brief Surface json decoding from already decoded bounds
//
// @tparam surface_t surface type
// @tparam bounds_t bounds base type the surface stores
//
// @param j json encoding of the surface
// @param bounds the bounds decoded from the surface's "bounds" payload
//
// @return shared pointer to the decoded surface
template <typename surface_t, typename bounds_t>
std::shared_ptr<Acts::Surface> surfaceFromJsonT(
    const nlohmann::json& j,
    std::shared_ptr<const Acts::SurfaceBounds> bounds) {
  Acts::Transform3 sTransform =
      Acts::Transform3JsonConverter::fromJson(j["transform"]);

  std::shared_ptr<const bounds_t> sBounds{};
  if (bounds != nullptr) {
    sBounds = std::dynamic_pointer_cast<const bounds_t>(bounds);
    // The surface/bounds pairing is only checked here, at decode time
    if (sBounds == nullptr) {
      throw std::invalid_argument(
          "Surface of type " + j.at("type").get<std::string>() +
          " cannot be built from bounds of type " + boundsKey(bounds->type()));
    }
  }
  return Acts::Surface::makeShared<surface_t>(sTransform, std::move(sBounds));
}

// @brief Surface json decoding for surfaces that carry no bounds
//
// @tparam surface_t surface type
//
// @param j json encoding of the surface
// @param bounds must be nullptr for these surfaces
//
// @return shared pointer to the decoded surface
template <typename surface_t>
std::shared_ptr<Acts::Surface> boundlessSurfaceFromJsonT(
    const nlohmann::json& j,
    std::shared_ptr<const Acts::SurfaceBounds> bounds) {
  if (bounds != nullptr) {
    throw std::invalid_argument("Surface of type " +
                                j.at("type").get<std::string>() +
                                " does not accept bounds");
  }
  Acts::Transform3 sTransform =
      Acts::Transform3JsonConverter::fromJson(j["transform"]);
  return Acts::Surface::makeShared<surface_t>(sTransform);
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
  cfg.surfaceEncoder.registerFunction(surfaceToJsonT<PointSurface>);

  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<EllipseBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<RectangleBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<TrapezoidBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<DiamondBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<AnnulusBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(surfaceBoundsToJsonT<RadialBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<DiscTrapezoidBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<CylinderBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(surfaceBoundsToJsonT<ConeBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(surfaceBoundsToJsonT<LineBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(surfaceBoundsToJsonT<PointBounds>);
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<InfiniteBounds>);
  // ConvexPolygonBounds is templated on the vertex count N. Registering the
  // abstract ConvexPolygonBoundsBase lets a single encoder handle every
  // specialization; the runtime dynamic_cast inside TypeDispatcher resolves
  // them all to the same handler.
  cfg.surfaceBoundsEncoder.registerFunction(
      surfaceBoundsToJsonT<ConvexPolygonBoundsBase>);

  // Bounds decoders, keyed on BoundsType alone
  cfg.surfaceBoundsDecoder
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eEllipse),
                    boundsFromJsonT<EllipseBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eRectangle),
                    boundsFromJsonT<RectangleBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eTrapezoid),
                    boundsFromJsonT<TrapezoidBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eDiamond),
                    boundsFromJsonT<DiamondBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eAnnulus),
                    boundsFromJsonT<AnnulusBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eDisc),
                    boundsFromJsonT<RadialBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eDiscTrapezoid),
                    boundsFromJsonT<DiscTrapezoidBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eCylinder),
                    boundsFromJsonT<CylinderBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eCone),
                    boundsFromJsonT<ConeBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::eLine),
                    boundsFromJsonT<LineBounds>)
      .registerKind(boundsKey(SurfaceBounds::BoundsType::ePoint),
                    boundsFromJsonT<PointBounds>);

  // Boundless surfaces carry no bounds object to build
  cfg.surfaceBoundsDecoder.registerKind(
      boundsKey(SurfaceBounds::BoundsType::eBoundless),
      [](const nlohmann::json& /*j*/) -> std::shared_ptr<const SurfaceBounds> {
        return nullptr;
      });

  // The vertex count of a ConvexPolygonBounds is recoverable only from the
  // length of the values array, so decode into a dynamically sized polygon.
  cfg.surfaceBoundsDecoder.registerKind(
      boundsKey(SurfaceBounds::BoundsType::eConvexPolygon),
      [](const nlohmann::json& j) -> std::shared_ptr<const SurfaceBounds> {
        std::vector<double> bVector = j["values"];
        if (bVector.size() < 6 || bVector.size() % 2 != 0) {
          throw std::invalid_argument(
              "Invalid ConvexPolygonBounds 'values' array: need an even "
              "number of entries (>= 6) encoding at least 3 (x, y) vertices");
        }
        std::vector<Vector2> vertices;
        vertices.reserve(bVector.size() / 2);
        for (std::size_t i = 0; i < bVector.size(); i += 2) {
          vertices.emplace_back(bVector[i], bVector[i + 1]);
        }
        return std::make_shared<const ConvexPolygonBounds<PolygonDynamic>>(
            vertices);
      });

  // Surface decoders, keyed on SurfaceType alone. Each narrows the decoded
  // bounds to the base type it stores, so any bounds deriving from that base
  // works without being enumerated here.
  cfg.surfaceDecoder
      .registerKind(surfaceKey(Surface::SurfaceType::Plane),
                    surfaceFromJsonT<PlaneSurface, PlanarBounds>)
      .registerKind(surfaceKey(Surface::SurfaceType::Disc),
                    surfaceFromJsonT<DiscSurface, DiscBounds>)
      .registerKind(surfaceKey(Surface::SurfaceType::Cylinder),
                    surfaceFromJsonT<CylinderSurface, CylinderBounds>)
      .registerKind(surfaceKey(Surface::SurfaceType::Cone),
                    surfaceFromJsonT<ConeSurface, ConeBounds>)
      .registerKind(surfaceKey(Surface::SurfaceType::Straw),
                    surfaceFromJsonT<StrawSurface, LineBounds>)
      .registerKind(surfaceKey(Surface::SurfaceType::Point),
                    surfaceFromJsonT<PointSurface, PointBounds>)
      .registerKind(surfaceKey(Surface::SurfaceType::Perigee),
                    boundlessSurfaceFromJsonT<PerigeeSurface>);

  return cfg;
}

Acts::SurfaceJsonConverter::Config Acts::SurfaceJsonConverter::m_cfg =
    Acts::SurfaceJsonConverter::Config::defaultConfig();

std::shared_ptr<Acts::Surface> Acts::SurfaceJsonConverter::fromJson(
    const nlohmann::json& j) {
  std::shared_ptr<const SurfaceBounds> bounds =
      m_cfg.surfaceBoundsDecoder(j.at("bounds"));
  std::shared_ptr<Acts::Surface> mutableSf =
      m_cfg.surfaceDecoder(j, std::move(bounds));

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
  jSurface["bounds"] = m_cfg.surfaceBoundsEncoder(surface.bounds());
  return jSurface;
}
