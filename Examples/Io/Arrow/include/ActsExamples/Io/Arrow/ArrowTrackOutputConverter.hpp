// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Parquet/ArrowOutputConverter.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <memory>
#include <string>
#include <vector>

#include <arrow/api.h>

namespace ActsExamples {

/// Convert a @c ConstTrackContainer to an @c arrow::Table.
///
/// The output table has one row per event with list-valued columns for the
/// perigee parameters (d0, z0, phi, theta, qop), the majority truth particle
/// id, the per-track measurement (hit) indices, and a running track index.
/// The @c ParquetWriter stamps the @c event_id column.
class ACTS_ARROW_EXPORT ArrowTrackOutputConverter final
    : public ArrowOutputConverter {
 public:
  struct Config {
    /// Input @c ConstTrackContainer on the whiteboard.
    std::string inputTracks;
    /// Optional input track-to-particle matching on the whiteboard. If empty,
    /// @c majority_particle_id is filled with zeros.
    std::string inputTrackParticleMatching;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "tracks";
  };

  explicit ArrowTrackOutputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const override;

 private:
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};

  WriteDataHandle<std::shared_ptr<arrow::Table>> m_outputTable{
      this, "OutputTable"};
};

}  // namespace ActsExamples
