// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::detail {
/// Optional internal navigation adapter; existing navigators need no new API.
template <typename Navigator>
struct NavigationStepGuard {
  template <typename... Args>
  static void beforeStep(Args&&...) {}
  template <typename... Args>
  static void afterStep(Args&&...) {}
  template <typename... Args>
  static void targetStatus(Args&&...) {}
};
}  // namespace Acts::detail
