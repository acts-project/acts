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
#include "ActsFatras/Kernel/detail/physics_list_implementation.hpp"
#include "ActsFatras/Kernel/detail/process_signature_check.hpp"

namespace ActsFatras {

/// This is the PhysicsList struct that is used for fast simulation
///
/// Users can add a variable list of processes in order to drive the
/// physics simulation
///
/// The dependency on generator, detector and particle are templated
template <typename... processes>
struct PhysicsList : private Acts::detail::Extendable<processes...> {
 private:
  using Acts::detail::Extendable<processes...>::tuple;

 public:
  using Acts::detail::Extendable<processes...>::get;

  /// Call operator that broadcasts the call to the tuple()
  /// members of the list
  ///
  /// @tparam generator_t is the random number generator type
  /// @tparam detector_t is the detector information type used
  /// @tparam particle_t is the particle type used in simulation
  ///
  /// @param[in] gen is the generator object
  /// @param[in] det is the necessary detector information
  /// @param[in] in is the ingoing particle (can be modified)
  /// @param[in,out] out are the (eventually) outgoing particles
  ///
  /// @return indicator which would trigger an abort
  template <typename generator_t, typename detector_t, typename particle_t>
  bool operator()(generator_t &gen, const detector_t &det, particle_t &in,
                  std::vector<particle_t> &out) const {
    // clang-format off
    static_assert(Acts::detail::all_of_v<detail::process_signature_check_v<processes, generator_t, detector_t, particle_t>...>,
                  "not all processes support the specified interface");
    // clang-format on

    // create an emtpy particle vector
    typedef detail::physics_list_impl<processes...> impl;
    return impl::process(tuple(), gen, det, in, out);
  }
};

}  // namespace ActsFatras
