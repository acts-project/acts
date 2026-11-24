// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <concepts>
#include <utility>

namespace Acts::detail {

/// Call a function when the scope is left, by return or by exception. For
/// restoring state that a function changes only temporarily.
///
/// @tparam function_t A callable without arguments that does not throw
template <std::invocable function_t>
class ScopeExit {
 public:
  /// Constructor
  /// @param function Called on destruction
  explicit ScopeExit(function_t function) : m_function(std::move(function)) {}

  ~ScopeExit() noexcept { m_function(); }

  ScopeExit(const ScopeExit&) = delete;
  ScopeExit& operator=(const ScopeExit&) = delete;
  ScopeExit(ScopeExit&&) = delete;
  ScopeExit& operator=(ScopeExit&&) = delete;

 private:
  function_t m_function;
};

}  // namespace Acts::detail
