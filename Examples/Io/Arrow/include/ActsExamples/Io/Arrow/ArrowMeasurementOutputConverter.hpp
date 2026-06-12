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
#include "ActsExamples/EventData/Cluster.hpp"
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
#include <vector>

namespace ActsExamples {

/// Convert measurements into the per-event RECO table (`tracker_hits`): one
/// entry per measurement, the entry's position being the measurement id
/// referenced by track `hit_ids`. Carries the measured local parameters and
/// variances (always expanded to the full bound layout, with a `subspace`
/// bitmask saying which components were actually measured: bit0 = loc0,
/// bit1 = loc1, bit2 = time), the reco global position, cluster-shape
/// features from the geometric digitization, and truth LINKS only:
/// `particle_ids` (rows in the particle table) and `simhit_ids` (rows in the
/// `tracker_simhits` table of the same event, see
/// `ArrowSimHitOutputConverter`).
class ACTS_ARROW_EXPORT ArrowMeasurementOutputConverter final
    : public ArrowOutputConverter {
 public:
  struct Config {
    /// Measurement container (required). One output row per measurement.
    std::string inputMeasurements;
    /// Optional cluster container parallel to the measurements (produced by
    /// the geometric digitization). When unset - e.g. smearing-only
    /// digitization - the shape columns are emitted as zeros so the schema
    /// stays stable.
    std::string inputClusters;
    /// Optional input @c SimHitContainer; must be set together with
    /// @c inputSimHitMeasurementsMap. When unset (e.g. data, or a reco-only
    /// conversion) the @c particle_ids / @c simhit_ids link columns are
    /// emitted as empty lists.
    std::string inputSimHits;
    /// Optional input particle container used to resolve a contributing
    /// sim-hit's particle barcode to a row index in the corresponding parquet
    /// table. Must be the same container the @c ArrowParticleOutputConverter
    /// consumes for that table — leaving it empty forces the unmatched
    /// sentinel.
    std::string inputParticles;
    /// Optional sim-hit → measurement(s) map keyed by @c SimHitIndex; must be
    /// set together with @c inputSimHits. It is inverted internally to
    /// measurement → contributing sim-hits.
    std::string inputSimHitMeasurementsMap;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "tracker_hits";
    /// Tracking geometry used to find surfaces by @c GeometryIdentifier when
    /// projecting each measurement's local parameters back to global
    /// coordinates (required).
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Diagnostic only: sim-hits contributing to no measurement are absent
    /// from this table's links (they remain in the truth table). A warning is
    /// logged when the fraction exceeds this value.
    double maxUnmatchedSimHitFraction = 0.001;
    /// Resolves the @c detector subsystem id for a given geometry id.
    /// Defaults to reading the geometry id's @c extra byte.
    std::function<std::uint8_t(Acts::GeometryIdentifier)> detectorResolver =
        [](Acts::GeometryIdentifier gid) {
          return static_cast<std::uint8_t>(gid.extra());
        };
  };

  explicit ArrowMeasurementOutputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const override;

 private:
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{
      this, "InputMeasurements"};
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<SimHitMeasurementsMap> m_inputSimHitMeasurementsMap{
      this, "InputSimHitMeasurementsMap"};

  WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_outputTable{
      this, "OutputTable"};
};

}  // namespace ActsExamples
