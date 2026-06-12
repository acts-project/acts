// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

/// Load a ColliderML geometry ID map from a Parquet file.
///
/// Expected columns (with header): detector, volume, layer, surface,
/// acts_geo_id. The @c acts_geo_id value is parsed as a hex or decimal
/// unsigned 64-bit integer (e.g. @c 0x1000000000010001).
///
/// @param path  Path to the Parquet file.
/// @return Map from packed ColliderML key to @c Acts::GeometryIdentifier.
std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>
loadColliderMLGeoIdMap(const std::filesystem::path& path);

/// Convert ColliderML Arrow tables to ACTS EDM types.
///
/// Reads two Arrow tables placed on the whiteboard by @c ParquetReader —
/// one for particles and one for tracker hits — and emits any combination
/// of @c SimParticleContainer, @c SimHitContainer, and
/// @c MeasurementContainer depending on which output keys are non-empty.
///
/// Modes (controlled by which output keys are set in the config):
///   - Particles only:  set @c outputParticles, leave hits keys empty.
///   - + SimHits:       also set @c outputSimHits.
///   - + Measurements:  also set @c outputMeasurements; requires
///                      @c trackingGeometry and @c digiConfig.
///
/// @note SimHit momentum fields are zero-filled; ColliderML does not
///       record per-hit momentum.
class ColliderMLInputConverter : public IAlgorithm {
 public:
  struct Config {
    /// Whiteboard key for the particles Arrow table (from ParquetReader).
    std::string inputParticlesTable;
    /// Whiteboard key for the tracker-hits Arrow table (from ParquetReader).
    std::string inputHitsTable;

    /// Output key for @c SimParticleContainer. Empty = skip.
    std::string outputParticles;
    /// Output key for @c SimHitContainer. Empty = skip.
    std::string outputSimHits;
    /// Output key for @c MeasurementContainer. Empty = skip.
    std::string outputMeasurements;
    /// Output key for @c MeasurementSubset covering all measurements (required
    /// by CKF / SpacePointMaker). Empty = skip.
    std::string outputMeasurementSubset;
    /// Output key for @c MeasurementSimHitsMap. Empty = skip.
    std::string outputMeasSimHitsMap;
    /// Output key for @c MeasurementParticlesMap (measurement→particle).
    /// Empty = skip.
    std::string outputMeasParticlesMap;
    /// Output key for @c ParticleMeasurementsMap (particle→measurements).
    /// Required by TruthTrackFinder / TruthSeeding. Empty = skip.
    std::string outputParticleMeasurementsMap;

    /// Required when @c outputMeasurements is non-empty.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// Digitisation config used to determine subspace and covariance per
    /// surface. Required when @c outputMeasurements is non-empty.
    /// Load with @c readDigiConfigFromJson (acts.examples.json module).
    /// The hierarchy map is queried per hit with full fallback
    /// (sensitive → layer → volume).
    DigiConfigContainer digiConfig;

    /// ColliderML (det, vol, layer, surf) → ACTS GeometryIdentifier.
    /// Optional. When empty, geometry IDs are constructed directly from the
    /// ColliderML (volume, layer, surface) fields — only correct when the data
    /// was produced from a geometry whose ID scheme matches the current build.
    /// Load with @c loadColliderMLGeoIdMap().
    std::unordered_map<std::uint64_t, Acts::GeometryIdentifier> geoIdMap;

    /// Euclidean boundary tolerance (mm) for projecting ColliderML 3D hit
    /// positions onto sensor surfaces. ColliderML hits are full 3D positions
    /// inside the sensor volume; at incidence angles the projected 2D
    /// coordinate can land a few mm outside the exact sensor boundary.
    /// Default 5 mm tolerates physical incidence effects while still catching
    /// wrong-surface assignments from a stale geoIdMap (tens of mm off).
    double hitBoundsTolerance = 5.0;
  };

  ColliderMLInputConverter(const Config& cfg,
                           std::unique_ptr<const Acts::Logger> logger);

  ColliderMLInputConverter(const Config& cfg, Acts::Logging::Level level);

  ~ColliderMLInputConverter() override;

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  /// Fallback map built at construction when geoIdMap is empty.
  /// Keys are GeometryIdentifiers with only volume/layer/sensitive set;
  /// values are the full GeometryIdentifier (including extra) from the
  /// tracking geometry.
  std::unordered_map<Acts::GeometryIdentifier, Acts::GeometryIdentifier>
      m_volLaySenMap;

  ReadDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_inputParticles{
      this, "InputParticles"};
  ReadDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_inputHits{this,
                                                                 "InputHits"};

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};
  WriteDataHandle<SimHitContainer> m_outputSimHits{this, "OutputSimHits"};
  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "OutputMeasurements"};
  WriteDataHandle<MeasurementSubset> m_outputMeasurementSubset{
      this, "OutputMeasurementSubset"};
  WriteDataHandle<MeasurementSimHitsMap> m_outputMeasSimHitsMap{
      this, "OutputMeasSimHitsMap"};
  WriteDataHandle<MeasurementParticlesMap> m_outputMeasParticlesMap{
      this, "OutputMeasParticlesMap"};
  WriteDataHandle<ParticleMeasurementsMap> m_outputParticleMeasurementsMap{
      this, "OutputParticleMeasurementsMap"};
};

}  // namespace ActsExamples
