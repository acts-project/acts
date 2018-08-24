// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
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
  /// @brief check and correct parameter values
  ///
  /// @tparam params template parameter pack containing the multiple parameter
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
  template <ParID_t... params>
  struct value_corrector;

  /// @cond
  template <typename R, ParID_t... params>
  struct value_corrector_impl;

  template <ParID_t... params>
  struct value_corrector
  {
    using ParVector_t = ActsVector<ParValue_t, sizeof...(params)>;

    static void
    result(ParVector_t& values)
    {
      value_corrector_impl<ParVector_t, params...>::calculate(values, 0);
    }
  };

  template <typename R, ParID_t first, ParID_t... others>
  struct value_corrector_impl<R, first, others...>
  {
    static void
    calculate(R& values, unsigned int pos)
    {
      using parameter_type = typename par_type<first>::type;
      if (parameter_type::may_modify_value) {
        values(pos) = parameter_type::getValue(values(pos));
      }
      value_corrector_impl<R, others...>::calculate(values, pos + 1);
    }
  };

  template <typename R, ParID_t last>
  struct value_corrector_impl<R, last>
  {
    static void
    calculate(R& values, unsigned int pos)
    {
      using parameter_type = typename par_type<last>::type;
      if (parameter_type::may_modify_value) {
        values(pos) = parameter_type::getValue(values(pos));
      }
    }
  };
  /// @endcond
}  // namespace detail
/// @endcond
}  // namespace Acts
