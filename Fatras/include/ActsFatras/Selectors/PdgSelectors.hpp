// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsFatras {

template <int pdg_t>
struct AbsPdgSelector {
  // absolute Pdg selection
  const int saPDG = pdg_t;

  /// Return true for all particles with | pdg | matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.pdg() * particle.pdg() == saPDG * saPDG);
  }
};

template <int pdg_t>
struct PdgSelector {
  // Pdg selection
  const int saPDG = pdg_t;

  /// Return true for all particles with pdg matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.pdg() == saPDG);
  }
};

template <int pdg_t>
struct AbsPdgExcluder {
  // absolute Pdg selection
  const int saPDG = pdg_t;

  /// Return true for all particles with | pdg | matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return !(particle.pdg() * particle.pdg() == saPDG * saPDG);
  }
};

template <int pdg_t>
struct PdgExcluder {
  // Pdg selection
  const int saPDG = pdg_t;

  /// Return true for all particles with pdg matching
  /// the selection criteria
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return !(particle.pdg() == saPDG);
  }
};

}  // namespace ActsFatras
