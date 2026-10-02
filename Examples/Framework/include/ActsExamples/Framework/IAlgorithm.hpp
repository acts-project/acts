// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// Event processing algorithm interface.
///
/// This class provides default implementations for most interface methods and
/// and adds a default logger that can be used directly in subclasses.
/// Algorithm implementations only need to implement the `execute` method.
class IAlgorithm : public SequenceElement {
 public:
  /// Construct an algorithm with a name and optional logger.
  /// @param name The algorithm name
  /// @param logger The logger for this algorithm
  explicit IAlgorithm(const std::string& name,
                      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// The algorithm name.
  /// @return The name supplied to the constructor.
  std::string name() const override;

  /// Execute the algorithm for one event.
  ///
  /// This function must be implemented by subclasses.
  /// @param context The current event context.
  /// @return The processing status for this event.
  virtual ProcessCode execute(const AlgorithmContext& context) const = 0;

  /// Internal execute method forwards to the algorithm execute method as const
  /// @param context The algorithm context
  /// @return The processing status from execute().
  ProcessCode internalExecute(const AlgorithmContext& context) final;

  /// Initialize the algorithm
  /// @return Success by default.
  ProcessCode initialize() override { return ProcessCode::SUCCESS; }
  /// Finalize the algorithm
  /// @return Success by default.
  ProcessCode finalize() override { return ProcessCode::SUCCESS; }

  /// Return the type name used in debug output.
  /// @return "Algorithm".
  std::string_view typeName() const override { return "Algorithm"; }

 protected:
  /// Return the algorithm logger.
  /// @return The logger owned by this algorithm.
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  std::string m_name;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
