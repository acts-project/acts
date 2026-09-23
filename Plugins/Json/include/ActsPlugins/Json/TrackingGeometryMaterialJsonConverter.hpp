// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Utilities/TypeDispatcher.hpp"
#include "ActsPlugins/Json/JsonKindDispatcher.hpp"

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace Acts {

/// Versioned material-document conversion, independent of the legacy material
/// converters. Applies no material to geometry. Essential validation is always
/// performed by the reader; no runtime JSON Schema validator is required.
/// Version 1 supports surface material only; volume assignments are rejected.
class TrackingGeometryMaterialJsonConverter {
 public:
  /// Shared grid slab allocation.
  using SlabStore = std::shared_ptr<std::vector<MaterialSlab>>;
  /// Per-document encoding state, also available to application encoders.
  class EncodeContext {
   public:
    /// Register a store, preserving sharing by allocation identity.
    /// @param store Non-null allocation
    /// @return Stable name within this document
    std::string storeId(const SlabStore& store);

   private:
    friend class TrackingGeometryMaterialJsonConverter;
    std::map<std::string, SlabStore, std::less<std::string>> m_stores;
  };
  /// Per-document decoding state, also available to application decoders.
  class DecodeContext {
   public:
    /// Resolve a shared slab store by its document name.
    /// @param name Store name in the material payload
    /// @return Shared allocation, preserving identity across payloads
    /// @throws std::invalid_argument if the store is unknown
    SlabStore store(const std::string& name) const;

   private:
    friend class TrackingGeometryMaterialJsonConverter;
    std::map<std::string, SlabStore, std::less<std::string>> m_stores;
  };
  /// Concrete surface encoder dispatch.
  using SurfaceEncoder =
      TypeDispatcher<ISurfaceMaterial, nlohmann::json(EncodeContext&)>;
  /// Surface decoder dispatch by kind.
  using SurfaceDecoder =
      JsonKindDispatcher<std::unique_ptr<const ISurfaceMaterial>,
                         const DecodeContext&>;
  /// Built-in or application-extended serialization registry.
  struct Config {
    /// Surface encoders.
    SurfaceEncoder encodeSurface;
    /// Surface decoders.
    SurfaceDecoder decodeSurface{"kind", "surface material"};
    /// Register all supported built-in kinds.
    /// @return Default registry, extensible by the caller
    static Config defaultConfig();
  };
  /// File output options. Version is a property of this codec, not an option.
  struct Options {
    /// Text indentation.
    unsigned indentation{4};
    /// Compression level, used only for zstd output.
    int compressionLevel{9};
    /// Default options for optional arguments.
    /// @return Default output options
    static Options defaultOptions() { return {}; }
  };

  /// Construct with default or application-extended dispatchers.
  /// @param config Payload conversion registry
  explicit TrackingGeometryMaterialJsonConverter(
      Config config = Config::defaultConfig());
  /// Encode a complete material document.
  /// @param material Container, including an optional description
  /// @throws std::invalid_argument if volume assignments are present
  /// @return Versioned JSON document
  nlohmann::json toJson(const TrackingGeometryMaterial& material) const;
  /// Decode and validate a complete material document.
  /// @param encoded Versioned document, never the legacy hierarchy-map layout
  /// @return Assignments and optional description, without applying them
  TrackingGeometryMaterial fromJson(const nlohmann::json& encoded) const;
  /// Write JSON/CBOR with optional zstd, selected by the filename extension.
  /// @param material Material container
  /// @param path Output filename
  /// @param options Output formatting/compression
  void toFile(const TrackingGeometryMaterial& material,
              const std::filesystem::path& path,
              const Options& options = Options::defaultOptions()) const;
  /// Read a document, detecting encoding/compression from file contents.
  /// @param path Input filename
  /// @return Decoded assignments, without applying them
  TrackingGeometryMaterial fromFile(const std::filesystem::path& path) const;

 private:
  Config m_config;
};

}  // namespace Acts
