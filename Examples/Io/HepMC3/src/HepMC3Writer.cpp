// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

ActsExamples::HepMC3WriterAscii::HepMC3WriterAscii(const Config&& cfg,
                                                   Acts::Logging::Level lvl)
    : WriterT(cfg.inputCollection, "HepMC3EventWriter", lvl), m_cfg(cfg) {
  if (m_cfg.outputStemFileName.empty())
    throw std::invalid_argument("Missing output stem file name");
}

ActsExamples::ProcessCode ActsExamples::HepMC3WriterAscii::writeT(
    const ActsExamples::AlgorithmContext& ctx,
    const std::vector<std::shared_ptr<HepMC3::GenEvent>>& events) {
  auto path = perEventFilepath(
      m_cfg.outputDir, m_cfg.outputStemFileName + ".hepmc3", ctx.eventNumber);
  std::cout << "Path: " << path << std::endl;
  HepMC3::WriterAscii writer(path);

  for (const auto& event : events)
    writer.write_event(*event);

  writer.close();
  return writer.failed() ? ActsExamples::ProcessCode::ABORT
                         : ActsExamples::ProcessCode::SUCCESS;
}
