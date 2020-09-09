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

template <typename R, BoundIndices... kIndices>
struct ValueCorrectorImpl;

template <typename R, BoundIndices kFirst, BoundIndices... kOthers>
struct ValueCorrectorImpl<R, kFirst, kOthers...> {
  static void correct(R& values, unsigned int pos) {
    using Traits = ParameterTraits<BoundIndices, kFirst>;

    if (Traits::may_modify_value) {
      values(pos) = Traits::getValue(values(pos));
    }
    ValueCorrectorImpl<R, kOthers...>::calculate(values, pos + 1u);
  }
};

template <typename R, BoundIndices kLast>
struct ValueCorrectorImpl<R, kLast> {
  static void correct(R& values, unsigned int pos) {
    using Traits = ParameterTraits<BoundIndices, kLast>;

    if (Traits::may_modify_value) {
      values(pos) = Traits::getValue(values(pos));
    }
  }
};

/// Correct parameter values.
///
/// @tparam kIndices template parameter pack containing the multiple parameter
///                identifiers
///
/// Values in the given vector are interpreted as values for the given
/// parameters. As those they are checked whether they are inside the allowed
/// range and corrected if necessary.
///
/// Invocation:
///   - value_corrector<params...>::result(parVector) where `parVector`
///     contains `sizeof...(params)` elements
///
/// @post All values in the argument `parVector` are within the valid
///       parameter range.
template <BoundIndices... kIndices>
struct ValueCorrector {
  using ParametersVector = ActsVector<BoundScalar, sizeof...(kIndices)>;

  static void correct(ParametersVector& values) {
    ValueCorrectorImpl<ParametersVector, kIndices...>::correct(values, 0u);
  }
};

}  // namespace detail
}  // namespace Acts
