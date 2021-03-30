// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/PrealignedDetectorElement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <mutex>
#include <vector>

namespace ActsExamples {

namespace Contextual {

/// @brief A mockup service that rotates the modules in a
/// simple tracking geometry
///
/// It acts on the PrealingedDetector element, i.e. all alignment
/// constants are already in memory and a random number is thrown
/// to pick one randomly for this event.
class PrealignedDecorator : public IContextDecorator {
 public:
  using LayerStore = std::vector<std::shared_ptr<PrealignedDetectorElement>>;
  using DetectorStore = std::vector<LayerStore>;

  /// @brief nested configuration struct
  struct Config {
    /// The detector store (filled at creation creation)
    DetectorStore detectorStore;

    /// The number of iovs
    unsigned int nIovs = 100;

    std::shared_ptr<RandomNumbers> randomNumberSvc = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param logger The logging framework
  PrealignedDecorator(const Config& cfg,
                      std::unique_ptr<const Acts::Logger> logger =
                          Acts::getDefaultLogger("PrealignedDecorator",
                                                 Acts::Logging::INFO));

  virtual ~PrealignedDecorator() = default;

  /// @brief decorates (adds, modifies) the AlgorithmContext
  /// with a new IOV context
  ///
  /// @param context the bare (or at least non-const) Event context
  ProcessCode decorate(AlgorithmContext& context) final override;

  /// @brief decorator name() for screen output
  const std::string& name() const final override { return m_name; }

 private:
  Config m_cfg;                                  ///< the configuration class
  std::unique_ptr<const Acts::Logger> m_logger;  ///!< the logging instance
  std::string m_name = "PrealignedDecorator";

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace Contextual
}  // namespace ActsExamples
