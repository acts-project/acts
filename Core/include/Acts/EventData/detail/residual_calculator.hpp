// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// Acts include(s)
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
/// @cond detail
namespace detail {
/// @brief calculate residuals from two parameter vectors
///
/// @tparam parameter_indices_t Specifying the underlying parameter index enum
/// @tparam params template parameter pack containing the multiple parameter
///                identifiers
///
/// Calculate the difference between the two given vectors with parameter
/// values. Possible corrections for bounded or cyclic parameters are applied.
///
/// @return residual_calculator<params...>::result(first,second) yields the
///         residuals of `first` with respect to `second`
template <typename parameter_indices_t, parameter_indices_t... params>
struct residual_calculator;
/// @cond

template <typename parameter_indices_t, typename R,
          parameter_indices_t... params>
struct residual_calculator_impl;

template <typename parameter_indices_t, parameter_indices_t... params>
struct residual_calculator {
  using ParVector_t = ActsVector<ParValue_t, sizeof...(params)>;

  static ParVector_t result(const ParVector_t& test, const ParVector_t& ref) {
    ParVector_t result;
    residual_calculator_impl<parameter_indices_t, ParVector_t,
                             params...>::calculate(result, test, ref, 0);
    return result;
  }
};

template <typename parameter_indices_t, typename R, parameter_indices_t first,
          parameter_indices_t... others>
struct residual_calculator_impl<parameter_indices_t, R, first, others...> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) =
        ParameterTypeFor<parameter_indices_t, first>::type::getDifference(
            test(pos), ref(pos));

    residual_calculator_impl<parameter_indices_t, R, others...>::calculate(
        result, test, ref, pos + 1);
  }
};

template <typename parameter_indices_t, typename R, parameter_indices_t last>
struct residual_calculator_impl<parameter_indices_t, R, last> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) =
        ParameterTypeFor<parameter_indices_t, last>::type::getDifference(
            test(pos), ref(pos));
  }
};
/// @endcond
}  // namespace detail
/// @endcond
}  // namespace Acts