// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
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

/// Convert simulated tracker hits into the per-event TRUTH table
/// (`tracker_simhits`): one entry per sim-hit, ALL sim-hits in container
/// order. The entry's position is the sim-hit id that the measurement table
/// (`ArrowMeasurementOutputConverter`) references via its `simhit_ids`
/// column. The table is standalone-complete for re-digitization: position +
/// time, 4-momentum at the hit, deposited energy, particle link, hit index
/// along the trajectory, and sensor identification.
class ACTS_ARROW_EXPORT ArrowSimHitOutputConverter final
    : public ArrowOutputConverter {
 public:
  struct Config {
    /// Input @c SimHitContainer on the whiteboard (required).
    std::string inputSimHits;
    /// Optional input particle container used to resolve a sim-hit's particle
    /// barcode to a row index in the corresponding parquet table. Must be the
    /// same container the @c ArrowParticleOutputConverter consumes for that
    /// table — leaving it empty forces the unmatched sentinel.
    std::string inputParticles;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "tracker_simhits";
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

  WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_outputTable{
      this, "OutputTable"};
};

}  // namespace ActsExamples
