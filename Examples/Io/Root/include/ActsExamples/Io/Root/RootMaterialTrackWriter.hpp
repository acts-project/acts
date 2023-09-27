// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <mutex>
#include <unordered_map>
class TFile;
class TTree;

namespace Acts {
// Using some short hands for Recorded Material
using RecordedMaterial = MaterialInteractor::result_type;
// And recorded material track
// - this is start:  position, start momentum
//   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;
}  // namespace Acts

namespace ActsExamples {

/// @class RootMaterialTrackWriter
///
/// @brief Writes out MaterialTrack collections from a root file
///
/// This service is the root implementation of the IWriterT.
/// It writes out a MaterialTrack which is usually generated from
/// Geant4 material mapping
class RootMaterialTrackWriter
    : public WriterT<std::unordered_map<size_t, Acts::RecordedMaterialTrack>> {
 public:
  struct Config {
    std::string collection =
        "material-tracks";                     ///< material collection to write
    std::string filePath = "";                 ///< path of the output file
    std::string fileMode = "RECREATE";         ///< file access mode
    std::string treeName = "material-tracks";  ///< name of the output tree

    /// Re-calculate total values from individual steps (for cross-checks)
    bool recalculateTotals = false;
    /// Write aut pre and post step (for G4), otherwise central step position
    bool prePostStep = false;
    /// Write the surface to which the material step correpond
    bool storeSurface = false;
    /// Write the volume to which the material step correpond
    bool storeVolume = false;
    /// Collapse consecutive interactions of a single surface into a single
    /// interaction
    bool collapseInteractions = false;
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level logging level
  RootMaterialTrackWriter(const Config& config,
                          Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootMaterialTrackWriter() override;

  /// Framework intialize method
  ActsExamples::ProcessCode finalize() override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 protected:
  // This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param clusters is the data to be written out
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::unordered_map<size_t, Acts::RecordedMaterialTrack>&
          materialtracks) override;

 private:
  /// The config class
  Config m_cfg;
  /// mutex used to protect multi-threaded writes
  std::mutex m_writeMutex;
  /// The output file name
  TFile* m_outputFile = nullptr;
  /// The output tree name
  TTree* m_outputTree = nullptr;

  /// Event identifier.
  uint32_t m_eventId = 0;

  float m_v_x = 0;    ///< start global x
  float m_v_y = 0;    ///< start global y
  float m_v_z = 0;    ///< start global z
  float m_v_px = 0;   ///< start global momentum x
  float m_v_py = 0;   ///< start global momentum y
  float m_v_pz = 0;   ///< start global momentum z
  float m_v_phi = 0;  ///< start phi direction
  float m_v_eta = 0;  ///< start eta direction
  float m_tX0 = 0;    ///< thickness in X0/L0
  float m_tL0 = 0;    ///< thickness in X0/L0

  std::vector<float> m_step_sx;      ///< step x (start) position (optional)
  std::vector<float> m_step_sy;      ///< step y (start) position (optional)
  std::vector<float> m_step_sz;      ///< step z (start) position (optional)
  std::vector<float> m_step_x;       ///< step x position
  std::vector<float> m_step_y;       ///< step y position
  std::vector<float> m_step_z;       ///< step z position
  std::vector<float> m_step_ex;      ///< step x (end) position (optional)
  std::vector<float> m_step_ey;      ///< step y (end) position (optional)
  std::vector<float> m_step_ez;      ///< step z (end) position (optional)
  std::vector<float> m_step_dx;      ///< step x direction
  std::vector<float> m_step_dy;      ///< step y direction
  std::vector<float> m_step_dz;      ///< step z direction
  std::vector<float> m_step_length;  ///< step length
  std::vector<float> m_step_X0;      ///< step material x0
  std::vector<float> m_step_L0;      ///< step material l0
  std::vector<float> m_step_A;       ///< step material A
  std::vector<float> m_step_Z;       ///< step material Z
  std::vector<float> m_step_rho;     ///< step material rho

  std::vector<std::uint64_t>
      m_sur_id;  ///< ID of the suface associated with the step
  std::vector<int32_t>
      m_sur_type;              ///< Type of the suface associated with the step
  std::vector<float> m_sur_x;  ///< x position of the center of the suface
                               ///< associated with the step
  std::vector<float> m_sur_y;  ///< y position of the center of the suface
                               ///< associated with the step
  std::vector<float> m_sur_z;  ///< z position of the center of the suface
                               ///< associated with the step
  std::vector<float>
      m_sur_pathCorrection;  ///< path correction when associating
                             ///< material to the given surface
  std::vector<float>
      m_sur_range_min;  ///< Min range of the suface associated with the step
  std::vector<float>
      m_sur_range_max;  ///< Max range of the suface associated with the step

  std::vector<std::uint64_t>
      m_vol_id;  ///< ID of the volume associated with the step
};

}  // namespace ActsExamples
