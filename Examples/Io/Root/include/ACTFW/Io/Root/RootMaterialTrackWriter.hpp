// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <mutex>

#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WriterT.hpp"

class TFile;
class TTree;

namespace Acts {
// Using some short hands for Recorded Material
using RecordedMaterial = MaterialInteractor::result_type;
// And recorded material track
// - this is start:  position, start momentum
//   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3D, Acts::Vector3D>, RecordedMaterial>;
}  // namespace Acts

namespace FW {

/// @class RootMaterialTrackWriter
///
/// @brief Writes out MaterialTrack collections from a root file
///
/// This service is the root implementation of the IWriterT.
/// It writes out a MaterialTrack which is usually generated from
/// Geant4 material mapping
class RootMaterialTrackWriter
    : public WriterT<std::vector<Acts::RecordedMaterialTrack>> {
 public:
  struct Config {
    std::string collection =
        "material-tracks";                     ///< material collection to write
    std::string filePath = "";                 ///< path of the output file
    std::string fileMode = "RECREATE";         ///< file access mode
    std::string treeName = "material-tracks";  ///< name of the output tree
    TFile* rootFile = nullptr;                 ///< common root file

    /// Re-calculate total values from individual steps (for cross-checks)
    bool recalculateTotals = false;
    /// Write aut pre and post step (for G4), otherwise central step position
    bool prePostStep = false;
    /// Write the surface to which the material step correpond
    bool storesurface = false;
  };

  /// Constructor with
  /// @param cfg configuration struct
  /// @param output logging level
  RootMaterialTrackWriter(const Config& cfg,
                          Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~RootMaterialTrackWriter() override;

  /// Framework intialize method
  FW::ProcessCode endRun() final override;

 protected:
  // This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ///
  /// @param ctx The Algorithm context with per event information
  /// @param clusters is the data to be written out
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::vector<Acts::RecordedMaterialTrack>&
                         materialtracks) final override;

 private:
  /// The config class
  Config m_cfg;
  /// mutex used to protect multi-threaded writes
  std::mutex m_writeMutex;
  /// The output file name
  TFile* m_outputFile;
  /// The output tree name
  TTree* m_outputTree;

  float m_v_x;    ///< start global x
  float m_v_y;    ///< start global y
  float m_v_z;    ///< start global z
  float m_v_px;   ///< start global momentum x
  float m_v_py;   ///< start global momentum y
  float m_v_pz;   ///< start global momentum z
  float m_v_phi;  ///< start phi direction
  float m_v_eta;  ///< start eta direction
  float m_tX0;    ///< thickness in X0/L0
  float m_tL0;    ///< thickness in X0/L0

  std::vector<float> m_step_sx;      ///< step x (start) position (optional)
  std::vector<float> m_step_sy;      ///< step y (start) position (optional)
  std::vector<float> m_step_sz;      ///< step z (start) position (optional)
  std::vector<float> m_step_x;       ///< step x position
  std::vector<float> m_step_y;       ///< step y position
  std::vector<float> m_step_z;       ///< step z position
  std::vector<float> m_step_ex;      ///< step x (end) position (optional)
  std::vector<float> m_step_ey;      ///< step y (end) position (optional)
  std::vector<float> m_step_ez;      ///< step z (end) position (optional)
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
      m_sur_range_min;  ///< Min range of the suface associated with the step
  std::vector<float>
      m_sur_range_max;  ///< Max range of the suface associated with the step
};

}  // namespace FW
