// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

#include <G4Track.hh>

namespace ActsExamples {

/// The G4SteppingAction that is called for every step in
/// the simulation process.
///
/// It checks whether a sensitive volume is present (via string tag)
/// and records (if necessary) the hit.
class G4TrackSelector {
 public:
  /// Configuration of the Stepping action
  struct Config {
    /// Selection for hit recording
    bool charged = true;
    bool neutral = false;
    bool primary = true;
    bool secondary = false;
  };

  G4TrackSelector(Config cfg) : m_cfg(cfg) {}

  bool operator()(const G4Track& track) const {
    G4PrimaryParticle* primaryParticle =
        track.GetDynamicParticle()->GetPrimaryParticle();

    G4double absCharge =
        std::abs(track.GetParticleDefinition()->GetPDGCharge());
    if (not m_cfg.charged and absCharge > 0.) {
      return false;
    }

    // Bail out if neutral & configured to do so
    if (not m_cfg.neutral and absCharge == 0.) {
      return false;
    }

    // Bail out if it is a primary & configured to be ignored
    if (not m_cfg.primary and primaryParticle != nullptr) {
      return false;
    }

    // Bail out if it is a secondary & configured to be ignored
    if (not m_cfg.secondary and primaryParticle == nullptr) {
      return false;
    }

    return true;
  }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
