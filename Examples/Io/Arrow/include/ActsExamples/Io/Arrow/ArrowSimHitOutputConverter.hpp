// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Parquet/ArrowOutputConverter.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ActsExamples {

/// Convert measurements (with their contributing sim-hit truth) to an
/// @c arrow::Table.
///
/// The output table has one row per event with list-valued columns, and one
/// entry per MEASUREMENT. Measurements are emitted in @c MeasurementContainer
/// iteration order, so the entry's position in the per-event lists equals its
/// measurement index — downstream tables (e.g. tracks) reference measurements
/// by that index, and there is no separate measurement-id column. Reco
/// @c x,y,z + geometry are one value per measurement; the contributing
/// sim-hits' truth (@c particle_id, @c true_x,true_y,true_z, @c time) are
/// nested lists (outer = per measurement, inner = per contributing sim-hit).
///
/// Sim-hits that contribute to no measurement are dropped (expected to be a
/// negligible fraction); the dropped fraction is checked at runtime against
/// @c maxUnmatchedSimHitFraction.
class ACTS_ARROW_EXPORT ArrowSimHitOutputConverter final
    : public ArrowOutputConverter {
 public:
  struct Config {
    /// Input @c SimHitContainer on the whiteboard (source of truth positions
    /// and particle barcodes for the nested per-measurement columns).
    std::string inputSimHits;
    /// Optional input particle container used to resolve a contributing
    /// sim-hit's particle barcode to a row index in the corresponding parquet
    /// table. Must be the same container the @c ArrowParticleOutputConverter
    /// consumes for that table — leaving it empty forces the unmatched sentinel.
    std::string inputParticles;
    /// Measurement container (required). One output row per measurement.
    std::string inputMeasurements;
    /// Sim-hit → measurement(s) map keyed by @c SimHitIndex (required). It is
    /// inverted internally to measurement → contributing sim-hits.
    std::string inputSimHitMeasurementsMap;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "simhits";
    /// Tracking geometry used to find surfaces by @c GeometryIdentifier when
    /// projecting each measurement's local parameters back to global
    /// coordinates (required).
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Sim-hits contributing to no measurement are dropped from the table.
    /// If the dropped fraction exceeds this, conversion throws (default 0.1%).
    double maxUnmatchedSimHitFraction = 0.001;
    /// Resolves the @c detector subsystem id for a given hit's geometry id.
    /// Defaults to reading the geometry id's @c extra byte; users can swap
    /// in any custom mapping (e.g. by volume or by surface lookup) when the
    /// geometry-construction side hasn't stamped @c extra yet.
    std::function<std::uint8_t(Acts::GeometryIdentifier)> detectorResolver =
        [](Acts::GeometryIdentifier gid) {
          return static_cast<std::uint8_t>(gid.extra());
        };
  };

  explicit ArrowSimHitOutputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Build a resolver from a volume-id -> detector-id lookup table.
  ///
  /// This returns a pure C++ callable so Python can configure the mapping
  /// once without paying a Python callback roundtrip for each hit.
  static std::function<std::uint8_t(Acts::GeometryIdentifier)>
  makeVolumeIdDetectorResolver(
      const std::unordered_map<std::uint32_t, std::uint8_t>& volumeToDetector,
      std::uint8_t defaultValue = 255);

  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const override;

 private:
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  Config m_cfg;

  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<SimHitMeasurementsMap> m_inputSimHitMeasurementsMap{
      this, "InputSimHitMeasurementsMap"};

  WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_outputTable{
      this, "OutputTable"};
};

}  // namespace ActsExamples
