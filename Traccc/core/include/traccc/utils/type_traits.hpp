// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::details {

/// Helper trait for detecting when a type is a non-const version of another
///
/// This comes into play multiple times to enable certain constructors
/// conditionally through SFINAE.
///
template <typename CTYPE, typename NCTYPE>
struct is_same_nc {
  static constexpr bool value = false;
};

template <typename TYPE>
struct is_same_nc<TYPE, TYPE> {
  static constexpr bool value = true;
};

template <typename TYPE>
struct is_same_nc<const TYPE, TYPE> {
  static constexpr bool value = true;
};

}  // namespace traccc::details
