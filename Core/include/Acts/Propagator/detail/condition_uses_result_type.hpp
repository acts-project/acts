// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <type_traits>
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  template <typename T, bool has_action = true>
  struct condition_uses_result_type_impl
  {
    static constexpr bool value = has_result_type_v<action_type_t<T>>;
  };

  template <typename T>
  struct condition_uses_result_type_impl<T, false> : std::false_type
  {
  };

  template <typename T>
  struct condition_uses_result_type
      : condition_uses_result_type_impl<T, has_action_type_v<T>>
  {
  };

}  // namespace detail

}  // namespace Acts