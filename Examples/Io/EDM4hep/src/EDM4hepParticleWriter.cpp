// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepParticleWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>

namespace ActsExamples {

EDM4hepParticleWriter::EDM4hepParticleWriter(
    const EDM4hepParticleWriter::Config& cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputParticles, "EDM4hepParticleWriter", lvl), m_cfg(cfg) {
  // TODO
}

ProcessCode EDM4hepParticleWriter::writeT(
    const AlgorithmContext& ctx, const SimParticleContainer& particles) {
  // TODO

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
