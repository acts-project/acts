// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>

class TFile;
class TTree;

namespace ActsExamples {

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
  ~RootParticleWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] particles are the particle to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SimParticleContainer& particles) override;

 private:
  Config m_cfg;
  std::mutex m_writeMutex;
  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;
  /// Event identifier.
  uint32_t m_eventId = 0;
  /// Event-unique particle identifier a.k.a barcode.
  std::vector<uint64_t> m_particleId;
  /// Particle type a.k.a. PDG particle number
  std::vector<int32_t> m_particleType;
  /// Production process type, i.e. what generated the particle.
  std::vector<uint32_t> m_process;
  /// Production position components in mm.
  std::vector<float> m_vx;
  std::vector<float> m_vy;
  std::vector<float> m_vz;
  // Production time in ns.
  std::vector<float> m_vt;
  /// Total momentum in GeV
  std::vector<float> m_p;
  /// Momentum components in GeV.
  std::vector<float> m_px;
  std::vector<float> m_py;
  std::vector<float> m_pz;
  /// Mass in GeV.
  std::vector<float> m_m;
  /// Charge in e.
  std::vector<float> m_q;
  // Derived kinematic quantities
  /// Direction pseudo-rapidity.
  std::vector<float> m_eta;
  /// Direction angle in the transverse plane.
  std::vector<float> m_phi;
  /// Transverse momentum in GeV.
  std::vector<float> m_pt;
  // Decoded particle identifier; see Barcode definition for details.
  std::vector<uint32_t> m_vertexPrimary;
  std::vector<uint32_t> m_vertexSecondary;
  std::vector<uint32_t> m_particle;
  std::vector<uint32_t> m_generation;
  std::vector<uint32_t> m_subParticle;
};

}  // namespace ActsExamples
