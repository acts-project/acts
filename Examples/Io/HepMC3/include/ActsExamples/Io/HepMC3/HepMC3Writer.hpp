// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"
#include <string>
#include <HepMC/GenEvent.h>
#include <HepMC3/WriterAscii.h>

namespace ActsExamples {

/// HepMC3 event writer.
struct HepMC3WriterAscii final
    : public WriterT<std::vector<std::shared_ptr<HepMC3::GenEvent>>> {
  struct Config {
    // The output directory
    std::string outputDir = "";
    // The stem of output file names
    std::string outputStemFileName = "";
    // The input collection
    std::string inputCollection = "";
  };

  /// @brief Constructor
  ///
  /// @param [in] cfg Config of the writer
  /// @param [in] lvl The level of the logger
  HepMC3WriterAscii(const Config&& cfg, Acts::Logging::Level lvl);

  /// @brief Writing method
  ///
  /// @param [in] ctx The context of this algorithm
  /// @param [in] events The recorded HepMC3 events
  ///
  /// @return Code describing whether the writing was successful
  ProcessCode writeT(const ActsExamples::AlgorithmContext& ctx,
                     const std::vector<std::shared_ptr<HepMC3::GenEvent>>&
                         events) final override;

 private:
  /// The configuration of this writer
  Config m_cfg;
};
}  // namespace ActsExamples
