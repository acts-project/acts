// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Parquet/ArrowOutputConverter.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <arrow/api.h>

namespace ActsExamples {

/// Convert a @c SimHitContainer to an @c arrow::Table.
///
/// The output table has one row per event with list-valued columns. Hits are
/// emitted in @c SimHitContainer iteration order, so the row index of a hit
/// inside the per-event list equals its @c SimHitIndex; downstream tables
/// (e.g. tracks) can therefore reference hits by that index.
///
/// When @c inputMeasurements, @c inputSimHitMeasurementsMap, and
/// @c trackingGeometry are all provided, the digitized cluster position
/// (local-to-global from the matched measurement) is written into @c x,y,z.
/// Otherwise those columns are filled with NaN. The truth position is always
/// written into @c true_x,true_y,true_z.
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
    /// Optional measurement container. Required (together with the two below)
    /// to fill the digitized @c x,y,z columns; otherwise those are NaN.
    std::string inputMeasurements;
    /// Optional sim-hit → measurement(s) inverse map; keyed by @c SimHitIndex.
    std::string inputSimHitMeasurementsMap;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "simhits";
    /// Tracking geometry used to find surfaces by @c GeometryIdentifier when
    /// projecting the matched measurement's local parameters back to global
    /// coordinates. Required for digitized @c x,y,z to be populated.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Map from @c GeometryIdentifier::volume() to a detector subsystem id
    /// (e.g. pixel=1, strip=2). Volumes absent from the map are written as 0.
    /// The original EDM4hep collection name is the natural source for this
    /// distinction but is lost in the internal EDM, so it has to be supplied
    /// explicitly.
    std::map<std::uint8_t, std::uint8_t> subsystemByVolume;
  };

  explicit ArrowSimHitOutputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

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

  WriteDataHandle<std::shared_ptr<arrow::Table>> m_outputTable{this,
                                                               "OutputTable"};
};

}  // namespace ActsExamples
