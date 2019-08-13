// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>
#include <utility>

#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
/// @cond
// forward declaration
template <ParID_t... params>
class ParameterSet;
template <FreeID_t... params>
class FreeParameterSet;
/// @endcond

/// @cond detail
namespace detail {
/// @brief generate ParameterSet type containing all defined parameters
///
/// @return full_parset<Policy>::type is equivalent to
///         `ParameterSet<Policy,ID_t(0),ID_t(1),...,ID_t(N-1)>` where @c ID_t
///         is a @c typedef to `Policy::par_id_type` and @c N is the total
///         number of parameters
struct full_parset {

  template <ParID_t v, typename C>
  struct add_to_value_container;

  template <ParID_t v, ParID_t... others>
  struct add_to_value_container<v, std::integer_sequence<ParID_t, others...>> {
    using type = std::integer_sequence<ParID_t, others..., v>;
  };

  template <typename T, unsigned int N>
  struct tparam_generator {
    using type = typename add_to_value_container<
        static_cast<T>(N), typename tparam_generator<T, N - 1>::type>::type;
  };

  template <typename T>
  struct tparam_generator<T, 0> {
    using type = std::integer_sequence<T, static_cast<T>(0)>;
  };

  template <typename T>
  struct converter;
  
  template <ParID_t... values>
  struct converter<std::integer_sequence<ParID_t, values...>> {
    using type = ParameterSet<values...>;
  };
  
  using type = typename converter<
      typename tparam_generator<ParID_t, BoundParsDim - 1>::type>::type;
};

struct full_freeparset {

  template <FreeID_t v, typename C>
  struct add_to_value_container;

  template <FreeID_t v, FreeID_t... others>
  struct add_to_value_container<v, std::integer_sequence<FreeID_t, others...>> {
    using type = std::integer_sequence<FreeID_t, others..., v>;
  };

  template <typename T, unsigned int N>
  struct tparam_generator {
    using type = typename add_to_value_container<
        static_cast<T>(N), typename tparam_generator<T, N - 1>::type>::type;
  };

  template <typename T>
  struct tparam_generator<T, 0> {
    using type = std::integer_sequence<T, static_cast<T>(0)>;
  };

  template <typename T>
  struct converter;

  template <FreeID_t... values>
  struct converter<std::integer_sequence<FreeID_t, values...>> {
    using type = FreeParameterSet<values...>;
  };
  
  using type = typename converter<
      typename tparam_generator<FreeID_t, FreeParsDim - 1>::type>::type;
};
}  // namespace detail
/// @endcond
}  // namespace Acts
