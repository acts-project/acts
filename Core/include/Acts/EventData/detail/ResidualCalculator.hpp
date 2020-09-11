// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/detail/ParameterTraits.hpp"

namespace Acts {
namespace detail {

template <typename indices_t, typename R, indices_t... kIndices>
struct ResidualCalculatorImpl;

template <typename indices_t, typename R, indices_t kFirst,
          indices_t... kOthers>
struct ResidualCalculatorImpl<indices_t, R, kFirst, kOthers...> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) =
        ParameterTraits<indices_t, kFirst>::getDifference(test(pos), ref(pos));
    ResidualCalculatorImpl<indices_t, R, kOthers...>::calculate(result, test,
                                                                ref, pos + 1u);
  }
};

template <typename indices_t, typename R, indices_t kLast>
struct ResidualCalculatorImpl<indices_t, R, kLast> {
  static void calculate(R& result, const R& test, const R& ref,
                        unsigned int pos) {
    result(pos) =
        ParameterTraits<indices_t, kLast>::getDifference(test(pos), ref(pos));
  }
};

/// @brief calculate residuals from two parameter vectors
///
/// @tparam indices_t Specifying the underlying parameter index enum
/// @tparam kIndices template parameter pack containing the multiple parameter
///                identifiers
///
/// Calculate the difference between the two given vectors with parameter
/// values. Possible corrections for bounded or cyclic parameters are applied.
///
/// @return ResidualCalculator<indices_t, kIndices...>::calculate(test, ref) yields the
///         residuals of `test` with respect to `ref`
template <typename indices_t, indices_t... kIndices>
struct ResidualCalculator {
  using ParametersVector = ActsVector<BoundScalar, sizeof...(kIndices)>;

  static ParametersVector calculate(const ParametersVector& test,
                                    const ParametersVector& ref) {
    ParametersVector residual;
    ResidualCalculatorImpl<indices_t, ParametersVector, kIndices...>::calculate(
        residual, test, ref, 0u);
    return residual;
  }
};

}  // namespace detail
}  // namespace Acts
