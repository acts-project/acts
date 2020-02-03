// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"
#include "ActsFatras/Kernel/detail/selector_list_implementation.hpp"
#include "ActsFatras/Kernel/detail/selector_signature_check.hpp"

namespace ActsFatras {

/// @brief This is the SelectorList struct that is used for fast simulation
///
/// Users can add a variable list of selectors in order to drive the
/// physics simulation. Selectors can access particle information and
/// detector information to decide whether a process is to take place
template <bool inclusive, typename... selectors>
struct SelectorListAXOR : private Acts::detail::Extendable<selectors...> {
 private:
  static_assert(not Acts::detail::has_duplicates_v<selectors...>,
                "same selector type specified several times");

  using Acts::detail::Extendable<selectors...>::tuple;

 public:
  using Acts::detail::Extendable<selectors...>::get;

  /// Call operator that is that broadcasts the call to the tuple()
  ///
  /// @tparam detector_t is the detector type used in simulation
  /// @tparam particle_t is the particle type used in simulation
  /// @tparam inclusive steers || (true) or && (false) combination
  ///
  /// @param[in] detector the current detector/material information
  /// @param[in] particle to be checked for further processing
  ///
  /// @return indicator if the particle is accepted
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &detector,
                  const particle_t &particle) const {
    // clang-format off
    static_assert(Acts::detail::all_of_v<detail::selector_list_signature_check_v<selectors, detector_t, particle_t>...>,
                  "not all particle selectors support the specified interface");
    // clang-format on

    // create an emtpy particle vector
    typedef detail::selector_list_impl<selectors...> impl;
    return impl::select(tuple(), detector, particle, inclusive);
  }
};

template <typename... selectors>
using SelectorListOR = SelectorListAXOR<true, selectors...>;

template <typename... selectors>
using SelectorListAND = SelectorListAXOR<false, selectors...>;

}  // namespace ActsFatras
