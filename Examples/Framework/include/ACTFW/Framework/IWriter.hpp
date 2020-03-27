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

#include <string>

#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// Event data writer interface.
///
/// Get data from the event store and write it to disk. The writer can have
/// internal state and implementations are responsible to handle concurrent
/// calls.
class IWriter {
 public:
  virtual ~IWriter() = default;

  /// The writer name.
  virtual std::string name() const = 0;

  /// Write data from one event.
  virtual ProcessCode write(const AlgorithmContext& context) = 0;

  /// End the run (e.g. aggregate statistics, write down output, close files).
  virtual ProcessCode endRun() = 0;
};

}  // namespace FW
