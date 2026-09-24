// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"
#include "ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp"

#include <filesystem>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

/// MaterialMapping output adapter for versioned JSON/CBOR surface material
/// maps, optionally compressed. Add it to
/// MaterialMapping::Config::materialWriters to serialize the finalized map
/// using Acts::TrackingGeometryMaterialJsonConverter.
///
/// writeMaterial() implements the mapping output interface. The separate
/// write() convenience method exports assignments directly from an existing
/// geometry.
class TrackingGeometryMaterialJsonWriter final : public IMaterialWriter {
 public:
  /// Writer configuration.
  struct Config {
    /// Full output path; the extension selects the encoding and compression.
    std::filesystem::path filePath;
    /// Add deferred one-bin grid placeholders for surfaces without material.
    bool includeNonMaterial = false;
    /// Precision, indentation and compression settings.
    Acts::TrackingGeometryMaterialJsonConverter::Options options;
  };

  /// Construct a writer.
  /// @param config Output configuration
  /// @param level Logging level
  TrackingGeometryMaterialJsonWriter(const Config& config,
                                     Acts::Logging::Level level);

  /// Write assignments. Volume material is unsupported and rejected.
  /// @param material Material assignments
  void writeMaterial(const Acts::TrackingGeometryMaterial& material) override;

  /// Export existing material assignments from geometry, preserving stable
  /// keys. Surfaces without material are omitted unless includeNonMaterial is
  /// set.
  /// @param gctx Geometry context used to resolve placeholder ranges
  /// @param geometry Geometry to export
  void write(const Acts::GeometryContext& gctx,
             const Acts::TrackingGeometry& geometry);

  /// @return Output configuration
  const Config& config() const { return m_config; }

 private:
  Config m_config;
};
}  // namespace ActsExamples
