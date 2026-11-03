// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

/// Command line options used to configure performance measurements
class performance : public interface {
 public:
  /// @name Options
  /// @{

  /// Whether to run performance checks
  bool run = false;

  /// @}

  /// Constructor
  performance();

  std::unique_ptr<configuration_printable> as_printable() const override;
};  // struct performance

}  // namespace traccc::opts
