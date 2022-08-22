// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <mutex>
#include <vector>

class TChain;

namespace ActsExamples {

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
    std::string treeName = "material-tracks";  ///< name of the output tree
    std::vector<std::string> fileList;         ///< List of input files

    /// Whether the events are ordered or not
    bool orderedEvents = true;
    // Read surface information for the root file
    bool readSurface = false;
  };

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  RootMaterialTrackReader(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootMaterialTrackReader();

  /// Framework name() method
  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(
      const ActsExamples::AlgorithmContext& context) final override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The logger
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  size_t m_events = 0;

  /// The batch size (number of track per events)
  size_t m_batch = 0;

  /// The input tree name
  TChain* m_inputChain = nullptr;

  /// Event identifier.
  uint32_t m_eventId;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

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

  std::vector<float>* m_step_x = new std::vector<float>;   ///< step x position
  std::vector<float>* m_step_y = new std::vector<float>;   ///< step y position
  std::vector<float>* m_step_z = new std::vector<float>;   ///< step z position
  std::vector<float>* m_step_dx = new std::vector<float>;  ///< step x direction
  std::vector<float>* m_step_dy = new std::vector<float>;  ///< step y direction
  std::vector<float>* m_step_dz = new std::vector<float>;  ///< step z direction
  std::vector<float>* m_step_length = new std::vector<float>;  ///< step length
  std::vector<float>* m_step_X0 = new std::vector<float>;  ///< step material x0
  std::vector<float>* m_step_L0 = new std::vector<float>;  ///< step material l0
  std::vector<float>* m_step_A = new std::vector<float>;   ///< step material A
  std::vector<float>* m_step_Z = new std::vector<float>;   ///< step material Z
  std::vector<float>* m_step_rho =
      new std::vector<float>;  ///< step material rho

  std::vector<std::uint64_t>* m_sur_id =
      new std::vector<std::uint64_t>;  ///< ID of the suface associated with the
                                       ///< step
  std::vector<float>* m_sur_x =
      new std::vector<float>;  ///< x position of the center of the suface
                               ///< associated with the step
  std::vector<float>* m_sur_y =
      new std::vector<float>;  ///< y position of the center of the suface
                               ///< associated with the step
  std::vector<float>* m_sur_z =
      new std::vector<float>;  ///< z position of the center of the suface
                               ///< associated with the step
};

}  // namespace ActsExamples
