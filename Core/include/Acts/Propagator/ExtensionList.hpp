// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/extension_list_implementation.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"

namespace Acts {

template <typename... extensions>
struct ExtensionList : private detail::Extendable<extensions...>
{
private:
  static_assert(not detail::has_duplicates_v<extensions...>,
                "same extension type specified several times");

  using detail::Extendable<extensions...>::tuple;

  using impl = detail::extension_list_impl<extensions...>;

public:
  using detail::Extendable<extensions...>::get;

  template <typename stepper_state_t>
  bool
  k(const stepper_state_t& state,
    Vector3D&              knew,
    const Vector3D&        bField,
    const int              i     = 0,
    const double           h     = 0,
    const Vector3D&        kprev = Vector3D())
  {
    return impl::k(tuple(), state, knew, bField, i, h, kprev);
  }

  template <typename stepper_state_t, typename stepper_data_t>
  bool
  finalize(stepper_state_t&      state,
           const double          h,
           const stepper_data_t& data,
           ActsMatrixD<7, 7>& D)
  {
    return impl::finalize(tuple(), state, h, data, D);
  }
};

}  // namespace Acts
