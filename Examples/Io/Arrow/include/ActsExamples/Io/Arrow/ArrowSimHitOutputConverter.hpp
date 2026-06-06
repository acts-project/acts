// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
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

/// Convert a @c SimHitContainer to an @c arrow::Table.
///
/// The output table has one row per event with list-valued columns. Hits are
/// emitted in @c SimHitContainer iteration order, so the row index of a hit
/// inside the per-event list equals its @c SimHitIndex; downstream tables
/// (e.g. tracks) can therefore reference hits by that index.
///
/// When @c inputClusters and @c inputSimHitMeasurementsMap are both provided,
/// the precomputed digitized cluster position (@c Cluster::globalPosition) of
/// the matched measurement is written into @c x,y,z. Clusters have a one-to-one
/// relation with measurements, so the @c SimHitMeasurementsMap (keyed by
/// @c SimHitIndex, valued by measurement index) doubles as a sim-hit → cluster
/// map. Otherwise those columns are filled with NaN. The truth position is
/// always written into @c true_x,true_y,true_z.
class ACTS_ARROW_EXPORT ArrowSimHitOutputConverter final
    : public ArrowOutputConverter {
 public:
  struct Config {
    /// Input @c SimHitContainer on the whiteboard.
    std::string inputSimHits;
    /// Optional input particle container used to resolve the hit's particle
    /// barcode to a row index in the corresponding parquet table. Must be the
    /// same container the @c ArrowParticleOutputConverter consumes for that
    /// table — leaving it empty forces the unmatched sentinel.
    std::string inputParticles;
    /// Optional cluster container. Required (together with the map below) to
    /// fill the digitized @c x,y,z columns from @c Cluster::globalPosition;
    /// otherwise those are NaN. Clusters are indexed one-to-one with
    /// measurements.
    std::string inputClusters;
    /// Optional sim-hit → measurement(s) inverse map; keyed by @c SimHitIndex.
    /// Because clusters and measurements share indices, the values double as
    /// cluster indices.
    std::string inputSimHitMeasurementsMap;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "simhits";
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
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};
  ReadDataHandle<SimHitMeasurementsMap> m_inputSimHitMeasurementsMap{
      this, "InputSimHitMeasurementsMap"};

  WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_outputTable{
      this, "OutputTable"};
};

}  // namespace ActsExamples
