// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvExaTrkXGraphWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <stdexcept>
#include <vector>

#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_namedtuple.hpp>

struct GraphData {
  int64_t edge0;
  int64_t edge1;
  float weight;
  DFE_NAMEDTUPLE(GraphData, edge0, edge1, weight);
};

ActsExamples::CsvExaTrkXGraphWriter::CsvExaTrkXGraphWriter(
    const ActsExamples::CsvExaTrkXGraphWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputGraph, "CsvExaTrkXGraphWriter", level),
      m_cfg(config) {}

ActsExamples::ProcessCode ActsExamples::CsvExaTrkXGraphWriter::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const std::pair<std::vector<int64_t>, std::vector<float>>& graph) {
  std::string path = perEventFilepath(
      m_cfg.outputDir, m_cfg.outputStem + ".csv", ctx.eventNumber);

  dfe::NamedTupleCsvWriter<GraphData> writer(path);

  const auto& [edges, weights] = graph;

  for (auto i = 0ul; i < weights.size(); ++i) {
    GraphData edge{};
    edge.edge0 = edges[2 * i];
    edge.edge1 = edges[2 * i + 1];
    edge.weight = weights[i];
    writer.append(edge);
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
