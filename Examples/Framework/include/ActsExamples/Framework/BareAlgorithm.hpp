// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-05-11 Initial version
/// @date 2017-07-27 Simplify interface
/// @author Andreas Salzburger
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>

namespace ActsExamples {

/// A helper class for users to implement framework algorithms
///
/// This class provides default implementations for most interface methods and
/// and adds a default logger that can be used directly in subclasses.
/// Algorithm implementations only need to implement the `execute` method.
class BareAlgorithm : public IAlgorithm {
 public:
  /// Constructor
  ///
  /// @name The algorithm name
  /// @level The logging level for this algorithm
  BareAlgorithm(std::string name,
                Acts::Logging::Level level = Acts::Logging::INFO);

  /// The algorithm name.
  std::string name() const final override;

  /// Execute the algorithm for one event.
  ///
  /// This function must be implemented by subclasses.
  virtual ProcessCode execute(
      const AlgorithmContext& context) const override = 0;

  /// Initialize the algorithm
  ProcessCode initialize() const override { return ProcessCode::SUCCESS; }
  /// Finalize the algorithm
  ProcessCode finalize() const override { return ProcessCode::SUCCESS; }

 protected:
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  std::string m_name;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
