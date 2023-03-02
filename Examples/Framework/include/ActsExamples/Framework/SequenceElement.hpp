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

/// Event processing interface.
///
class SequenceElement {
 public:
  virtual ~SequenceElement() = default;

  /// The algorithm name.
  virtual std::string name() const = 0;

  /// Initialize the algorithm
  virtual ProcessCode initialize() = 0;

  /// Finalize the algorithm
  virtual ProcessCode finalize() = 0;

  /// Internal method to execute the algorithm for one event.
  /// @note Usually, you should not override this method
  virtual ProcessCode internalExecute(const AlgorithmContext& context) = 0;
};

}  // namespace ActsExamples
