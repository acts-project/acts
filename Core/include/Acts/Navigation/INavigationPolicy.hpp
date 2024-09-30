// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Utilities/DelegateChainBuilder.hpp"

#include <type_traits>

namespace Acts {

class TrackingVolume;
class INavigationPolicy;

template <typename T>
concept NavigationPolicyConcept = requires {
  requires std::is_base_of_v<INavigationPolicy, T>;
  // Has a conforming update method
  requires requires(T policy, const NavigationArguments& args) {
    { policy.updateState(args) };
  };
};

class INavigationPolicy {
 public:
  static void noopUpdate(const NavigationArguments& /*unused*/) {}

  virtual ~INavigationPolicy() = default;

  virtual void connect(NavigationDelegate& delegate) const = 0;

 protected:
  template <NavigationPolicyConcept T>
  void connectDefault(NavigationDelegate& delegate) const {
    // This cannot be a concept because we use it in CRTP below
    const auto* self = static_cast<const T*>(this);
    DelegateChainBuilder{delegate}.add<&T::updateState>(self).store(delegate);
  }
};

}  // namespace Acts
