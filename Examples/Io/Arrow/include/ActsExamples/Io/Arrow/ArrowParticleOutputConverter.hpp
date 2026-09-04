// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Parquet/ArrowOutputConverter.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

/// Convert a @c SimParticleContainer to an @c arrow::Table.
///
/// The output table has one row per particle with columns for id, PDG code,
/// charge, mass, and the initial-state four-momentum / four-position. The
/// table is placed on the whiteboard under the configured key; the
/// @c ParquetWriter picks it up from there and stamps the @c event_id column.
class ACTS_ARROW_EXPORT ArrowParticleOutputConverter final
    : public ArrowOutputConverter {
 public:
  struct Config {
    /// Input @c SimParticleContainer on the whiteboard.
    std::string inputParticles;
    /// Output whiteboard key for the resulting @c arrow::Table.
    std::string outputTable = "particles";
  };

  explicit ArrowParticleOutputConverter(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~ArrowParticleOutputConverter() override;

  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const override;

 private:
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_outputTable{
      this, "OutputTable"};
};

}  // namespace ActsExamples
