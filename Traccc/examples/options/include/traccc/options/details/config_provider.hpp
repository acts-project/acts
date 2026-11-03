// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::opts {
/**
 * @brief Mixin type to indicate that some set of program options can be
 * converted to some configuration type.
 *
 * @tparam Config The config type to which this can be converted
 */
template <typename Config>
class config_provider {
 public:
  using config_type = Config;

  virtual operator config_type() const = 0;
};
}  // namespace traccc::opts
