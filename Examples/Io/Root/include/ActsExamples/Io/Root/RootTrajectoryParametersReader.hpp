// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
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

/// @class RootTrajectoryParametersReader
///
/// @brief Reads in TrackParameter information from a root file
/// and fills it into a Acts::BoundTrackParameter format
class RootTrajectoryParametersReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::string outputTracks = "outputTracks";  ///< track collection to read
    std::string outputParticles =
        "outputParticles";                        ///< track collection to read
    std::string treeName = "trackparams_fitter";  ///< name of the output tree
    std::string inputFile;  ///< The name of the input file
    std::string inputDir;   ///< The name of the input dir

    /// The default logger
    std::shared_ptr<const Acts::Logger> logger;

    /// The name of the service
    std::string name;

    /// Constructor
    /// @param lname The name of the reader
    /// @parqam lvl The log level for the logger
    Config(const std::string& lname = "TrackParameterReader",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  /// @param cfg The Configuration struct
  RootTrajectoryParametersReader(const Config& cfg);

  /// Destructor
  ~RootTrajectoryParametersReader();

  /// Framework name() method
  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(
      const ActsExamples::AlgorithmContext& context) final override;

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

  unsigned int m_eventNr{0};  ///< the event number
  std::vector<unsigned int>* m_multiTrajNr =
      new std::vector<unsigned int>;  ///< the multi-trajectory number
  std::vector<unsigned int>* m_subTrajNr =
      new std::vector<unsigned int>;  ///< the multi-trajectory sub-trajectory
                                      ///< number

  std::vector<unsigned long>* m_t_barcode =
      new std::vector<long unsigned int>;  ///< Truth particle barcode
  std::vector<int>* m_t_charge =
      new std::vector<int>;  ///< Truth particle charge
  std::vector<float>* m_t_time =
      new std::vector<float>;  ///< Truth particle time
  std::vector<float>* m_t_vx =
      new std::vector<float>;  ///< Truth particle vertex x
  std::vector<float>* m_t_vy =
      new std::vector<float>;  ///< Truth particle vertex y
  std::vector<float>* m_t_vz =
      new std::vector<float>;  ///< Truth particle vertex z
  std::vector<float>* m_t_px =
      new std::vector<float>;  ///< Truth particle initial momentum px
  std::vector<float>* m_t_py =
      new std::vector<float>;  ///< Truth particle initial momentum py
  std::vector<float>* m_t_pz =
      new std::vector<float>;  ///< Truth particle initial momentum pz
  std::vector<float>* m_t_theta =
      new std::vector<float>;  ///< Truth particle initial momentum theta
  std::vector<float>* m_t_phi =
      new std::vector<float>;  ///< Truth particle initial momentum phi
  std::vector<float>* m_t_pT =
      new std::vector<float>;  ///< Truth particle initial momentum pT
  std::vector<float>* m_t_eta =
      new std::vector<float>;  ///< Truth particle initial momentum eta

  std::vector<bool>* m_hasFittedParams =
      new std::vector<bool>;  ///< if the track has fitted parameter
  std::vector<float>* m_eLOC0_fit =
      new std::vector<float>;  ///< fitted parameter eBoundLoc0
  std::vector<float>* m_eLOC1_fit =
      new std::vector<float>;  ///< fitted parameter eBoundLoc1
  std::vector<float>* m_ePHI_fit =
      new std::vector<float>;  ///< fitted parameter ePHI
  std::vector<float>* m_eTHETA_fit =
      new std::vector<float>;  ///< fitted parameter eTHETA
  std::vector<float>* m_eQOP_fit =
      new std::vector<float>;  ///< fitted parameter eQOP
  std::vector<float>* m_eT_fit =
      new std::vector<float>;  ///< fitted parameter eT
  std::vector<float>* m_err_eLOC0_fit =
      new std::vector<float>;  ///< fitted parameter eLOC err
  std::vector<float>* m_err_eLOC1_fit =
      new std::vector<float>;  ///< fitted parameter eBoundLoc1 err
  std::vector<float>* m_err_ePHI_fit =
      new std::vector<float>;  ///< fitted parameter ePHI err
  std::vector<float>* m_err_eTHETA_fit =
      new std::vector<float>;  ///< fitted parameter eTHETA err
  std::vector<float>* m_err_eQOP_fit =
      new std::vector<float>;  ///< fitted parameter eQOP err
  std::vector<float>* m_err_eT_fit =
      new std::vector<float>;  ///< fitted parameter eT err
};

}  // namespace ActsExamples
