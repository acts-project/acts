// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
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
/// @tparam params template parameter pack containing the multiple parameter
///                identifiers
///
/// Calculate the difference between the two given vectors with parameter
/// values. Possible corrections for bounded or cyclic parameters are applied.
///
/// @return residual_calculator<params...>::result(first,second) yields the
///         residuals of `first` with respect to `second`
template <ParID_t... params>
struct residual_calculator;

template <FreeParametersIndices... params>
struct free_residual_calculator;

/// @cond
template <typename R, ParID_t... params>
struct residual_calculator_impl;

template <typename R, FreeParametersIndices... params>
struct free_residual_calculator_impl;

template <ParID_t... params>
struct residual_calculator {
  using ParVector_t = ActsVector<ParValue_t, sizeof...(params)>;

  static ParVector_t result(const ParVector_t& test, const ParVector_t& ref) {
    ParVector_t result;
    residual_calculator_impl<ParVector_t, params...>::calculate(result, test,
                                                                ref, 0);
    return result;
  }
};

template <FreeParametersIndices... params>
struct free_residual_calculator {
  using ParVector_t = ActsVector<ParValue_t, sizeof...(params)>;

  static ParVector_t result(const ParVector_t& test, const ParVector_t& ref) {
    ParVector_t result;
    free_residual_calculator_impl<ParVector_t, params...>::calculate(result, test,
                                                                ref, 0);
    return result;
  }
};

template <typename R, ParID_t first, ParID_t... others>
struct residual_calculator_impl<R, first, others...> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) = BoundParameterType<first>::getDifference(test(pos), ref(pos));
    residual_calculator_impl<R, others...>::calculate(result, test, ref,
                                                      pos + 1);
  }
};

template <typename R, FreeParametersIndices first, FreeParametersIndices... others>
struct free_residual_calculator_impl<R, first, others...> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) = FreeParameterType<first>::getDifference(test(pos), ref(pos));
    free_residual_calculator_impl<R, others...>::calculate(result, test, ref,
                                                      pos + 1);
  }
};

template <typename R, ParID_t last>
struct residual_calculator_impl<R, last> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) = BoundParameterType<last>::getDifference(test(pos), ref(pos));
  }
};

template <typename R, FreeParametersIndices last>
struct free_residual_calculator_impl<R, last> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) = FreeParameterType<last>::getDifference(test(pos), ref(pos));
  }
};
/// @endcond
}  // namespace detail
/// @endcond
}  // namespace Acts
