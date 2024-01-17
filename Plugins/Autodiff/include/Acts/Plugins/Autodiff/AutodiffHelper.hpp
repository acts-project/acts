// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace autodiff {
namespace detail {

template <typename To, typename From>
To copysign(To to, From &&from) {
  return (from >= 0 ? 1.0f : -1.0f) * to;
}

}  // namespace detail
}  // namespace autodiff
