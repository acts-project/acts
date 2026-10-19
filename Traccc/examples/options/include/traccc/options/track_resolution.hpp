// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/ambiguity_resolution/ambiguity_resolution_config.hpp"
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

/// Configuration for track ambiguity resolution
class track_resolution : public interface,
                         public config_provider<ambiguity_resolution_config> {
 public:
  /// Constructor
  track_resolution();

  /// Configuration conversion operators
  operator ambiguity_resolution_config() const override;

  /// Configuration conversion operators
  std::unique_ptr<configuration_printable> as_printable() const override;

 private:
  /// The internal configuration
  ambiguity_resolution_config m_config;
};  // class track_resolution

}  // namespace traccc::opts
