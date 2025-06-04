// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

#include <array>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace ActsExamples {

/// @brief Decorator for the alignment context
///
/// This decorator adds an alignment context to the algorithm context, allowing
/// algorithms to access alignment information.
///
/// This AlignmentDectorator can run in thwo modes:
///
/// 1. Read alignment data from an I/O store, which is a shared pointer to an
///    IAlignmentStore. The store is expected to provide contextual transforms
///    for surfaces.
/// 2. Generate alignment data on-the-fly using generators, which are functions
///    that take a reference to an `Acts::Transform3` and modify it. The
///    generators are expected to be provided with validity intervals (IOVs),
///    which are simple arrays with two elements representing the start and end
class AlignmentDecorator : public IContextDecorator {
 public:
  /// The interval of validity (IOV) for alignment data
  /// This is a simple array with two elements, representing the start and end
  /// of the validity interval.
  using IOV = std::array<std::size_t, 2>;

  /// Configuration struct
  struct Config {
    /// I/O mode - this if for showcase examples where the alignment data is
    /// read
    std::vector<std::tuple<IOV, std::shared_ptr<IAlignmentStore>>>
        iovStores;  //!< Stores for alignment data, each with a validity
                    //!< interval

    /// Generation mode - this is for showcase examples where the alignment data
    /// is generated`
    std::shared_ptr<IAlignmentStore> nominalStore =
        nullptr;  //!< The nominal alignment store

    /// Run garbage collection on the alignment data
    bool garbageCollection =
        false;  //!< Whether to run garbage collection on the alignment data

    /// Number of events to wait before garbage collection
    std::size_t gcInterval =
        100u;  //!< The number of events to wait before garbage collection

    std::vector<std::tuple<IOV, std::function<void(Acts::Transform3*)>>>
        iovGenerators;  //!< Generators for alignment data, each with a validity
                        //!< interval
  };

  /// Constructor with configuration and logger
  /// @param config Configuration for the alignment decorator
  /// @param level Logging level for the decorator
  AlignmentDecorator(const Config& config, Acts::Logging::Level level);

  /// Decorate the context with the alignment context
  /// @param context The algorithm context to decorate
  ProcessCode decorate(AlgorithmContext& context) override;

  /// Get the name of the decorator
  const std::string& name() const override { return m_name; }

 private:
  Config m_cfg;  //!< The configuration strcut

  std::unique_ptr<const Acts::Logger> m_logger;  ///!< The logging instance

  std::string m_name = "AlignmentDecorator";  //!< The name of the decorator

  std::mutex m_iovMutex;  //@< Mutex for IOV operations
  std::vector<std::tuple<IOV, std::shared_ptr<IAlignmentStore>, std::size_t>>
      m_iovStores;  //!< Stores for alignment data, each with a validity
                    //!< interval and a counter of events passed since last in
                    //!< use

  std::vector<std::tuple<IOV, std::function<void(Acts::Transform3*)>>>
      m_iovGenerators;  //!< Generators for alignment data, each with a validity
                        //!< interval

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
// namespace ActsExamples
