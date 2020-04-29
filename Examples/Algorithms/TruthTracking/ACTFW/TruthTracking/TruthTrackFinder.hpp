// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Framework/BareAlgorithm.hpp"

namespace FW {

/// Convert true particle tracks into "reconstructed" proto tracks.
///
/// For numbering consistency, this creates a proto track for each input
/// particle. Depending on the input particle selection it can contain zero
/// hits. This algorithm should be able to replace any other real track finder
/// in the reconstruction chain e.g. to validate algorithms further down
/// the chain.
class TruthTrackFinder final : public BareAlgorithm {
 public:
  struct Config {
    /// The input truth particles that should be used to create proto tracks.
    std::string inputParticles;
    /// The input hit-particles map collection.
    std::string inputHitParticlesMap;
    /// The output proto tracks collection.
    std::string outputProtoTracks;
  };

  TruthTrackFinder(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const override final;

 private:
  Config m_cfg;
};

}  // namespace FW
