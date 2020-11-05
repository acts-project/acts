// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include <vector>
#include <cstdint>
#include <string>

class TH1;

/// @brief This class transforms and stores a not necessarily normalised probability distribution into a cumulative probability distribution. The class allows the sampling of random numbers according the underlying distribution of the given probability distribution.
class TFCS1DFunctionInt32Histogram
{
  public:
	typedef uint32_t HistoContent_t;
    static constexpr HistoContent_t s_MaxValue = UINT32_MAX;
    
	/// @brief Constructor
	///
	/// @param [in] hist The histogram
	/// @param [in] name Name of this histogram
	/// @param [in] lvl The verbosity level of the logger
    TFCS1DFunctionInt32Histogram(const TH1* hist, std::string name, Acts::Logging::Level lvl);

    ///Function gets random number rnd in the range [0,1) as argument 
    ///and returns function value according to a histogram distribution
    double rndToFunction(double rnd) const;

    const std::vector<float>& getHistoBorders() const {return m_HistoBorders;};
    std::vector<float>& getHistoBorders() {return m_HistoBorders;};
    const std::vector<HistoContent_t>& getHistoContents() const {return m_HistoContents;};
    std::vector<HistoContent_t>& getHistoContents() {return m_HistoContents;};
      
  private:
  	/// @brief This function transforms a given distribution into a cumulative distribution
    void Initialize(const TH1* hist);
  
  const Acts::Logger& logger() const { return *m_logger; }
    
  /// Bin borders of the histogram
    std::vector<float> m_HistoBorders;
    /// Bin content of the cumulative distribution
    std::vector<HistoContent_t> m_HistoContents;
    /// The logger
    std::unique_ptr<const Acts::Logger> m_logger;
};