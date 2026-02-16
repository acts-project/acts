// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvSpacePointWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Io/Csv/CsvInputOutput.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <string>

#include "CsvOutputData.hpp"

namespace ActsExamples {

CsvSpacePointWriter::CsvSpacePointWriter(const Config& config,
                                         Acts::Logging::Level level)
    : WriterT(config.inputSpacepoints, "CsvSpacePointWriter", level),
      m_cfg(config) {}

CsvSpacePointWriter::~CsvSpacePointWriter() = default;

ProcessCode CsvSpacePointWriter::finalize() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ProcessCode CsvSpacePointWriter::writeT(
    const AlgorithmContext& ctx, const SimSpacePointContainer& spacepoints) {
  // Open per-event file for all components
  std::string pathSP =
      perEventFilepath(m_cfg.outputDir, "spacepoint.csv", ctx.eventNumber);

  BoostDescribeCsvWriter<SpacepointData> writerSP(pathSP,
                                                  m_cfg.outputPrecision);

  SpacepointData spData{};
  for (const auto& sp : spacepoints) {
    const auto slink1 = sp.sourceLinks()[0].get<IndexSourceLink>();
    spData.measurement_id_1 = slink1.index();
    spData.geometry_id_1 = slink1.geometryId().value();
    if (sp.sourceLinks().size() == 2) {
      const auto slink2 = sp.sourceLinks()[1].get<IndexSourceLink>();
      spData.measurement_id_2 = slink2.index();
      spData.geometry_id_2 = slink2.geometryId().value();
    } else {
      spData.measurement_id_2 = std::numeric_limits<std::uint64_t>::max();
      spData.geometry_id_2 = 0;  // invalid geoid
    }
    spData.x = sp.x() / Acts::UnitConstants::mm;
    spData.y = sp.y() / Acts::UnitConstants::mm;
    spData.z = sp.z() / Acts::UnitConstants::mm;
    spData.t = sp.t() ? *sp.t() / Acts::UnitConstants::ns
                      : std::numeric_limits<double>::quiet_NaN();
    spData.var_r = sp.varianceR() / Acts::UnitConstants::mm;
    spData.var_z = sp.varianceZ() / Acts::UnitConstants::mm;
    writerSP.append(spData);
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
