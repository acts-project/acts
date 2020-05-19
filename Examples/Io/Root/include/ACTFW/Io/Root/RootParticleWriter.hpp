// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <mutex>
#include <string>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WriterT.hpp"

class TFile;
class TTree;

namespace FW {

/// Write out particles as a flat TTree.
///
/// Each entry in the TTree corresponds to one particle for optimum writing
/// speed. The event number is part of the written data.
///
/// Safe to use from multiple writer threads. To avoid thread-saftey issues,
/// the writer must be the sole owner of the underlying file. Thus, the
/// output file pointer can not be given from the outside.
class RootParticleWriter final : public WriterT<SimParticleContainer> {
 public:
  struct Config {
    /// Input particle collection to write.
    std::string inputParticles;
    /// Path to the output file.
    std::string filePath;
    /// Output file access mode.
    std::string fileMode = "RECREATE";
    /// Name of the tree within the output file.
    std::string treeName = "particles";
  };

  /// Construct the particle writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  RootParticleWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Ensure underlying file is closed.
  ~RootParticleWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] particles are the particle to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SimParticleContainer& particles) final override;

 private:
  Config m_cfg;
  std::mutex m_writeMutex;
  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;
  /// Event identifier.
  uint32_t m_eventId;
  /// Event-unique particle identifier a.k.a barcode.
  uint64_t m_particleId;
  /// Particle type a.k.a. PDG particle number
  int32_t m_particleType;
  /// Production process type, i.e. what generated the particle.
  uint32_t m_process;
  /// Production position components in mm.
  float m_vx, m_vy, m_vz;
  // Production time in ns.
  float m_vt;
  /// Momentum components in GeV.
  float m_px, m_py, m_pz;
  /// Mass in GeV.
  float m_m;
  /// Charge in e.
  float m_q;
  // Derived kinematic quantities
  /// Direction pseudo-rapidity.
  float m_eta;
  /// Direction angle in the transverse plane.
  float m_phi;
  /// Transverse momentum in GeV.
  float m_pt;
  // Decoded particle identifier; see Barcode definition for details.
  uint32_t m_vertexPrimary;
  uint32_t m_vertexSecondary;
  uint32_t m_particle;
  uint32_t m_generation;
  uint32_t m_subParticle;
};

}  // namespace FW
