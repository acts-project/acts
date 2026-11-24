// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>
#include <utility>

#include <pybind11/pybind11.h>

/// Property-binding helpers shared by the track and track state proxy
/// bindings (EventDataTracks.cpp, EventDataTrackStates.cpp).
namespace ActsPython::EventDataHelpers {

/// The plain Eigen type that `access(self)` evaluates to, e.g.
/// `Eigen::Map<const Acts::BoundVector>::PlainObject` is `Acts::BoundVector`.
template <typename Proxy, typename Access>
using EigenValueT = typename std::decay_t<decltype(std::declval<Access>()(
    std::declval<const Proxy&>()))>::PlainObject;

/// Bind a property whose accessor returns an Eigen::Map, always read-only
/// (for accessors with no mutable overload, e.g. TrackStateProxy's "best
/// available" parameters()/covariance()).
template <typename Cls, typename Access>
void bindEigenReadonly(Cls& cls, const char* name, Access access) {
  using Proxy = typename Cls::type;
  using Value = EigenValueT<Proxy, Access>;
  cls.def_property_readonly(
      name, [access](const Proxy& self) { return Value{access(self)}; });
}

/// Bind a property whose accessor returns an Eigen::Map, read-write when
/// `!Proxy::ReadOnly`. `access` must be a generic lambda so the constness of
/// `self` selects the const/mutable overload.
template <typename Cls, typename Access>
void bindEigen(Cls& cls, const char* name, Access access) {
  using Proxy = typename Cls::type;
  using Value = EigenValueT<Proxy, Access>;
  auto get = [access](const Proxy& self) { return Value{access(self)}; };
  if constexpr (Proxy::ReadOnly) {
    cls.def_property_readonly(name, get);
  } else {
    cls.def_property(
        name, get, [access](Proxy& self, const Value& v) { access(self) = v; });
  }
}

/// Bind a property for an accessor overloaded on const (by value) vs. mutable
/// (by reference) with the same name, e.g. tipIndex, nMeasurements, chi2.
/// Read-write when `!Proxy::ReadOnly`; see @ref bindEigen for `access`.
template <typename Cls, typename Access>
void bindScalar(Cls& cls, const char* name, Access access) {
  using Proxy = typename Cls::type;
  using T = std::decay_t<decltype(access(std::declval<const Proxy&>()))>;
  auto get = [access](const Proxy& self) -> T { return access(self); };
  if constexpr (Proxy::ReadOnly) {
    cls.def_property_readonly(name, get);
  } else {
    cls.def_property(name, get,
                     [access](Proxy& self, T v) { access(self) = v; });
  }
}

}  // namespace ActsPython::EventDataHelpers
