// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2017-07-25
/// @author Moritz Kiehnn <msmk@cern.ch>

#pragma once

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"

#include <string>

namespace ActsExamples {

/// Event data writer interface.
///
/// Get data from the event store and write it to disk. The writer can have
/// internal state and implementations are responsible to handle concurrent
/// calls.
class IWriter : public SequenceElement {
 public:
  /// Write data from one event.
  virtual ProcessCode write(const AlgorithmContext& context) = 0;

  /// Internal execute method forwards to the write method as mutable
  /// @param context The algorithm context
  ProcessCode internalExecute(const AlgorithmContext& context) final {
    return write(context);
  }

  /// Fullfil the algorithm interface
  ProcessCode initialize() override { return ProcessCode::SUCCESS; }
};

}  // namespace ActsExamples
