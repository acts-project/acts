// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <mutex>
#include <vector>

#include "ACTFW/Framework/IReader.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

class TChain;

namespace FW {

/// @class RootMaterialTrackReader
///
/// @brief Reads in MaterialTrack information from a root file
/// and fills it into a format to be understood by the MaterialMapping
/// algorithm
class RootMaterialTrackReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::string collection =
        "material-tracks";                     ///< material collection to read
    std::string filePath = "";                 ///< path of the output file
    std::string treeName = "material-tracks";  ///< name of the output tree
    std::vector<std::string> fileList;         ///< The name of the input file

    unsigned int batchSize = 1;  ///!< The number of tracks per event

    /// The default logger
    std::shared_ptr<const Acts::Logger> logger;

    /// The name of the service
    std::string name;

    /// Constructor
    /// @param lname The name of the Material reader
    /// @parqam lvl The log level for the logger
    Config(const std::string& lname = "MaterialReader",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  /// @param cfg The Configuration struct
  RootMaterialTrackReader(const Config& cfg);

  /// Destructor
  ~RootMaterialTrackReader();

  /// Framework name() method
  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(const FW::AlgorithmContext& context) final override;

 private:
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }

  /// The config class
  Config m_cfg;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  size_t m_events = 0;

  /// The input tree name
  TChain* m_inputChain = nullptr;

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

  std::vector<float>* m_step_x = new std::vector<float>;  ///< step x position
  std::vector<float>* m_step_y = new std::vector<float>;  ///< step y position
  std::vector<float>* m_step_z = new std::vector<float>;  ///< step z position
  std::vector<float>* m_step_length = new std::vector<float>;  ///< step length
  std::vector<float>* m_step_X0 = new std::vector<float>;  ///< step material x0
  std::vector<float>* m_step_L0 = new std::vector<float>;  ///< step material l0
  std::vector<float>* m_step_A = new std::vector<float>;   ///< step material A
  std::vector<float>* m_step_Z = new std::vector<float>;   ///< step material Z
  std::vector<float>* m_step_rho =
      new std::vector<float>;  ///< step material rho
};

}  // namespace FW
