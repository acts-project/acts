// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

template <typename, typename H = void, DelegateType = DelegateType::NonOwning>
class DelegateChain;

/// This class is a thin wrapper around a vector of delegates that are supposed
/// to be called in order.
template <typename R, typename H, DelegateType O, typename... Args>
  requires(std::is_same_v<R, void>)  // Currently does not support return values
class DelegateChain<R(Args...), H, O> {
  using return_type = R;
  using holder_type = H;

 public:
  static constexpr DelegateType kOwnership = O;

  using DelegateType = Delegate<R(Args...), H, kOwnership>;

  /// Insert a delegate at the end of the chain
  /// @param delegate The delegate to append
  void push_back(DelegateType&& delegate) {
    assert(delegate.connected() && "Delegate is not connected");
    m_delegates.push_back(delegate);
  }

  /// Creates a new delegate at the end of the chain and immediately connects it
  /// with the given arguments.
  /// @tparam Callable The callable to connect
  /// @tparam Ts The argument types
  /// @param args The arguments to connect the delegate with
  template <auto Callable, typename... Ts>
  void connect_back(Ts&&... args) {
    auto& delegate = m_delegates.emplace_back();
    delegate.template connect<Callable>(std::forward<Ts>(args)...);
  }

  /// The number of delegates in the chain
  /// @return The number of delegates
  std::size_t size() const { return m_delegates.size(); }

  /// Check if the chain is empty
  /// @return True if the chain is empty, else false
  bool empty() const { return m_delegates.empty(); }

  /// The call operator that exposes the functionality of the @c DelegateChain type.
  /// @param args The arguments to call the contained functions with
  template <typename... Ts>
  void operator()(Ts&&... args) const {
    for (const auto& delegate : m_delegates) {
      delegate(std::forward<Ts>(args)...);
    }
  }

 private:
  std::vector<DelegateType> m_delegates;
};

}  // namespace Acts
