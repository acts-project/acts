// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

/// Event processing algorithm interface.
///
/// An algorithm must have no internal state and can communicate to the
/// rest of the world only by reading and writting to the event store.
class IAlgorithm {
 public:
  virtual ~IAlgorithm() = default;

  /// The algorithm name.
  virtual std::string name() const = 0;

  /// Execute the algorithm for one event.
  virtual ProcessCode execute(const AlgorithmContext& context) const = 0;
  /// Initialize the algorithm
  virtual ProcessCode initialize() const = 0;
  /// Finalize the algorithm
  virtual ProcessCode finalize() const = 0;
};

}  // namespace ActsExamples
