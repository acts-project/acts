// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

/// Option(s) for accelerator usage
class accelerator : public interface {
 public:
  /// @name Options
  /// @{

  /// Whether GPU texture memory should be used
  bool use_gpu_texture_memory = false;

  /// @}

  /// Constructor
  accelerator();

  std::unique_ptr<configuration_printable> as_printable() const override;
};  // struct accelerator

}  // namespace traccc::opts
