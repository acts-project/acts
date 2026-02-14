// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvVertexWriter.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>

#include "CsvOutputData.hpp"

namespace ActsExamples {

CsvVertexWriter::CsvVertexWriter(const Config& cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputVertices, "CsvVertexWriter", lvl), m_cfg(cfg) {
  // inputVertices is already checked by base constructor
  if (m_cfg.outputStem.empty()) {
    throw std::invalid_argument("Missing output filename stem");
  }
}

ProcessCode CsvVertexWriter::writeT(const AlgorithmContext& ctx,
                                    const VertexContainer& vertices) {
  auto pathVertices = perEventFilepath(
      m_cfg.outputDir, m_cfg.outputStem + ".csv", ctx.eventNumber);
  NamedTupleCsvWriter<VertexData> writer(pathVertices, m_cfg.outputPrecision);

  // Iterate over the vertices, and write out the 4 positions
  VertexData data;
  for (const auto& vertex : vertices) {
    data.x = vertex.fullPosition()[Acts::CoordinateIndices::eX];
    data.y = vertex.fullPosition()[Acts::CoordinateIndices::eY];
    data.z = vertex.fullPosition()[Acts::CoordinateIndices::eZ];
    data.T = vertex.fullPosition()[Acts::CoordinateIndices::eTime];

    writer.append(data);
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
