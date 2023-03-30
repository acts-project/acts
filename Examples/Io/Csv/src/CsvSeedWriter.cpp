// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSeedWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <ios>
#include <optional>
#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvSeedWriter::CsvSeedWriter(
    const ActsExamples::CsvSeedWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputSeeds, "CsvSeedWriter", level), m_cfg(config) {}

ActsExamples::CsvSeedWriter::~CsvSeedWriter() = default;

ActsExamples::ProcessCode ActsExamples::CsvSeedWriter::endRun() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::CsvSeedWriter::writeT(
    const AlgorithmContext& ctx, const SimSeedContainer& seeds) {
  // Open per-event file for all components
  std::string pathSP =
      perEventFilepath(m_cfg.outputDir, "seed.csv", ctx.eventNumber);

  dfe::NamedTupleCsvWriter<SeedData> writerSP(pathSP, m_cfg.outputPrecision);

  SeedData seedData{};
  for (const auto& seed : seeds) {
    const auto& spacepoints = seed.sp();

    const auto slink_1 =
        spacepoints[0]->sourceLinks()[0].get<IndexSourceLink>();
    const auto slink_2 =
        spacepoints[1]->sourceLinks()[0].get<IndexSourceLink>();
    const auto slink_3 =
        spacepoints[2]->sourceLinks()[0].get<IndexSourceLink>();

    seedData.measurement_id_1 = slink_1.index();
    seedData.geometry_id_1 = slink_1.geometryId().value();
    seedData.x_1 = spacepoints[0]->x();
    seedData.y_1 = spacepoints[0]->y();
    seedData.z_1 = spacepoints[0]->z();
    seedData.var_r_1 = spacepoints[0]->varianceR();
    seedData.var_z_1 = spacepoints[0]->varianceZ();

    seedData.measurement_id_2 = slink_2.index();
    seedData.geometry_id_2 = slink_2.geometryId().value();
    seedData.x_2 = spacepoints[1]->x();
    seedData.y_2 = spacepoints[1]->y();
    seedData.z_2 = spacepoints[1]->z();
    seedData.var_r_2 = spacepoints[1]->varianceR();
    seedData.var_z_2 = spacepoints[1]->varianceZ();

    seedData.measurement_id_3 = slink_3.index();
    seedData.geometry_id_3 = slink_3.geometryId().value();
    seedData.x_3 = spacepoints[2]->x();
    seedData.y_3 = spacepoints[2]->y();
    seedData.z_3 = spacepoints[2]->z();
    seedData.var_r_3 = spacepoints[2]->varianceR();
    seedData.var_z_3 = spacepoints[2]->varianceZ();

    seedData.z_vertex = seed.z();
    seedData.seed_quality = seed.seedQuality();

    writerSP.append(seedData);
  }
  return ActsExamples::ProcessCode::SUCCESS;
}
