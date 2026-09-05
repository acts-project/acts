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

namespace detail {
class MaterialJsonEncodeContext;
class MaterialJsonDecodeContext;
}  // namespace detail

/// @addtogroup json_plugin
/// @{

/// Static class performing the JSON conversion of surface material
///
/// The encoding side is a @c TypeDispatcher registered on the concrete
/// (or, for the grid family, the abstract templated) material types, the
/// decoding side a @c JsonKindDispatcher keyed on the payload type tag.
class SurfaceMaterialJsonConverter {
 public:
  /// Context collecting the slab stores of the document being written
  using EncodeContext = detail::MaterialJsonEncodeContext;

  /// Context carrying the slab stores of the document being read
  using DecodeContext = detail::MaterialJsonDecodeContext;

  /// Encoder type for the surface material
  using Encoder =
      TypeDispatcher<ISurfaceMaterial, nlohmann::json(EncodeContext&)>;

  /// Decoder type for the surface material
  using Decoder = JsonKindDispatcher<std::unique_ptr<const ISurfaceMaterial>,
                                     const DecodeContext&>;

  /// Configuration struct
  struct Config {
    /// Encoder for the surface material
    Encoder encoder{};

    /// Decoder for the surface material, keyed on the payload type tag
    Decoder decoder{jsonKey().typekey, "surface material"};

    /// Access the shared default configuration
    ///
    /// @return the default configuration instance
    static const Config& defaultConfig();
  };

  /// Delete the default constructor as the class is purely static
  SurfaceMaterialJsonConverter() = delete;

  /// Convert surface material into its json payload
  ///
  /// @param material the material to be converted
  /// @param config the converter configuration
  /// @param context the document context collecting the slab stores, or a
  ///        nullptr to inline any slab store and keep the payload
  ///        self-contained
  ///
  /// @return the json payload of the material, i.e. the value that goes
  ///         under the @c material key of a surface
  static nlohmann::json toJson(const ISurfaceMaterial& material,
                               const Config& config = Config::defaultConfig(),
                               EncodeContext* context = nullptr);

  /// Convert a json payload back into surface material
  ///
  /// @param jMaterial the json payload of the material
  /// @param config the converter configuration
  /// @param context the document context holding the slab stores, or a
  ///        nullptr if the payload is expected to be self-contained
  ///
  /// @return the decoded material, or a nullptr if the payload is flagged
  ///         as not participating in the material mapping
  static std::unique_ptr<const ISurfaceMaterial> fromJson(
      const nlohmann::json& jMaterial,
      const Config& config = Config::defaultConfig(),
      const DecodeContext* context = nullptr);
};

/// @}

}  // namespace Acts
