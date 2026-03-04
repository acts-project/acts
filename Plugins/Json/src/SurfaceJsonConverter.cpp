// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsPlugins/Json/DetrayJsonHelper.hpp"
#include "ActsPlugins/Json/SurfaceBoundsJsonConverter.hpp"

#include <memory>

#include <nlohmann/json_fwd.hpp>

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

  return cfg;
}

Acts::SurfaceJsonConverter::Config Acts::SurfaceJsonConverter::m_cfg =
    Acts::SurfaceJsonConverter::Config::defaultConfig();

std::shared_ptr<Acts::Surface> Acts::SurfaceJsonConverter::fromJson(
    const nlohmann::json& j) {
  return m_cfg.surfaceDecoder(j);
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
