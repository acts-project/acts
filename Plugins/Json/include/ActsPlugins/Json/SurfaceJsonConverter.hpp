// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/TypeDispatcher.hpp"
#include "ActsPlugins/Json/AlgebraJsonConverter.hpp"
#include "ActsPlugins/Json/JsonKindDispatcher.hpp"

#include <memory>
#include <tuple>

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

/// Static class performing JSON conversion of the surfaces
class SurfaceJsonConverter {
 public:
  /// Options for surface conversion
  struct Options {
    /// Transform serialization options
    Transform3JsonConverter::Options transformOptions =
        Transform3JsonConverter::Options{};
    /// Write material information
    bool writeMaterial = true;
    /// Write surface as portal
    bool portal = false;
  };

  /// Encoder type for the surface bounds
  using SurfaceBoundsEncoder = TypeDispatcher<SurfaceBounds, nlohmann::json()>;
  /// Encoder type for the surfaces
  using SurfaceEncoder =
      TypeDispatcher<Surface,
                     nlohmann::json(const GeometryContext&, const Options&)>;

  /// Deccoder type for the surfaces
  using SurfaceDecoder = JsonKindDispatcher<std::shared_ptr<Surface>>;

  /// Configuration struct
  struct Config {
    /// Encoder for the surfaces
    SurfaceEncoder surfaceEncoder{};
    /// Encoder for the surface bounds
    SurfaceBoundsEncoder surfaceBoundsEncoder{};

    /// Decoder for the surfaces
    SurfaceDecoder surfaceDecoder{};

    /// Default configuration construction
    ///
    /// @return default configuration
    static Config defaultConfig();
  };

  /// Delete the default constructor
  /// as the class is purely static (for now)
  SurfaceJsonConverter() = delete;

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

  /// @brief Set externally constructed configuration
  ///
  /// @note May be removed when the dispatcher migration is finished
  ///
  /// @param cfg configuration to use
  static void setConfig(const Config& cfg) { m_cfg = cfg; };

 private:
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

}  // namespace Acts
