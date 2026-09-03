// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Utilities/TypeDispatcher.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/GeometryJsonKeys.hpp"
#include "ActsPlugins/Json/JsonKindDispatcher.hpp"

#include <memory>

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// Static class performing the JSON conversion of surface material
///
/// The encoding side is a @c TypeDispatcher registered on the concrete
/// (or, for the grid family, the abstract templated) material types, the
/// decoding side a @c JsonKindDispatcher keyed on the payload type tag.
class SurfaceMaterialJsonConverter {
 public:
  /// Encoder type for the surface material
  using Encoder = TypeDispatcher<ISurfaceMaterial, nlohmann::json()>;

  /// Decoder type for the surface material
  using Decoder = JsonKindDispatcher<std::unique_ptr<const ISurfaceMaterial>>;

  /// Configuration struct
  struct Config {
    /// Encoder for the surface material
    Encoder encoder{};

    /// Decoder for the surface material, keyed on the payload type tag
    Decoder decoder{jsonKey().typekey, "surface material"};

    /// Default configuration construction
    ///
    /// @return default configuration
    static Config defaultConfig();
  };

  /// Delete the default constructor as the class is purely static
  SurfaceMaterialJsonConverter() = delete;

  /// Access the shared default configuration
  ///
  /// @return the default configuration instance
  static const Config& defaultConfig();

  /// Convert surface material into its json payload
  ///
  /// @param material the material to be converted
  /// @param config the converter configuration
  ///
  /// @return the json payload of the material, i.e. the value that goes
  ///         under the @c material key of a surface
  static nlohmann::json toJson(const ISurfaceMaterial& material,
                               const Config& config = defaultConfig());

  /// Convert a json payload back into surface material
  ///
  /// @param jMaterial the json payload of the material
  /// @param config the converter configuration
  ///
  /// @return the decoded material, or a nullptr if the payload is flagged
  ///         as not participating in the material mapping
  static std::unique_ptr<const ISurfaceMaterial> fromJson(
      const nlohmann::json& jMaterial, const Config& config = defaultConfig());
};

/// Convert surface material into the @c material entry of a json object
///
/// @param j Destination JSON object
/// @param material Source material, may be a nullptr
void to_json(nlohmann::json& j,
             const std::shared_ptr<const ISurfaceMaterial>& material);

/// Read surface material from the @c material entry of a json object
///
/// @param j Source JSON object
/// @param material Destination material
void from_json(const nlohmann::json& j,
               std::shared_ptr<const ISurfaceMaterial>& material);

/// @}

}  // namespace Acts
