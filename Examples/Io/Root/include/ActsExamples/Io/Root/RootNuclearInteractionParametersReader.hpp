// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <mutex>
#include <vector>

class TH1F;

namespace ActsExamples {

/// @brief This class reads data that is relevant for simulating nuclear interaction. Beside the plain reading, the data is transformed and stored for application in ActsFatras
class RootNuclearInteractionParametersReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::vector<std::string> fileList;         ///< The name of the input file
    std::string outputParametrisation = "parameters"; ///< The name of the stored parametrisation

    /// The default logger
    std::shared_ptr<const Acts::Logger> logger;

    /// The name of the service
    std::string name;

	/// Number of simulated events for normalising the nuclear interaction probability
	unsigned int nSimulatedEvents = 100;
	
    /// Constructor
    /// @param lname The name of the Material reader
    /// @param lvl The log level for the logger
    Config(const std::string& lname = "NuclearInteractionParameters",
           Acts::Logging::Level lvl = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  /// @param cfg The Configuration struct
  RootNuclearInteractionParametersReader(const Config& cfg);

  /// Destructor
  ~RootNuclearInteractionParametersReader();

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
 /// @brief This method builds decomposed cumulative probability distributions
/// out of a vector of proability distributions
///
/// @param [in] histos Vector of probability distributions
///
/// @return Vector containing the decomposed cumulative probability
/// distributions
std::vector<std::pair<std::vector<float>, std::vector<uint32_t>>> buildMaps(
    const std::vector<TH1F*>& histos) const;
    
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_cfg.logger; }

  /// The config class
  Config m_cfg;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;
};

}  // namespace ActsExamples
