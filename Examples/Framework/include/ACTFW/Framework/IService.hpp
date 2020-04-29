// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-05-23
/// @author Andreas Salburger
/// @author Moritz Kiehnn <msmk@cern.ch>

#pragma once

#include <string>

#include "ACTFW/Framework/AlgorithmContext.hpp"

namespace FW {

class WhiteBoard;

/// Service interface.
///
/// A service should be used to provide constant or slowly changing
/// per-event information, e.g. geometry with or without alignment, magnetic
/// field, ..., and to handle once-per-run tasks. In contrast to an
/// algorithm (i.e. inheriting from IAlgorithm), a service can have an internal
/// state and each
/// implementation has to ensure that concurrent calls are valid.
class IService {
 public:
  virtual ~IService() = default;

  /// The service name.
  virtual std::string name() const = 0;

  /// Start-of-run hook to be called before any events are processed.
  ///
  /// Should throw an exception for non-recoverable errors.
  virtual void startRun() = 0;

  /// Prepare per-event information.
  ///
  /// This is intended to add already existing information, e.g. geometry
  /// or conditions data, to the event store. While possible, complex
  /// operations should be better implemented as an regular algorithm.
  ///
  /// Should throw an exception on non-recoverable errors.
  virtual void prepare(AlgorithmContext& ctx) = 0;
};

}  // namespace FW
