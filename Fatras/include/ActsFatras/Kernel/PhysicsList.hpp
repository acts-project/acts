// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <bitset>
#include <tuple>

#include "Acts/Material/MaterialProperties.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// Combined set of physics processes and interactions for the simulation.
///
/// The physics processes are extendable by the user to be able to accomodate
/// the specific requirements. While the set of available physics processes must
/// be configured at compile-time, within that set processes can be selectively
/// disabled at run-time. By default all processes are applied.
template <typename... processes_t>
class PhysicsList {
 public:
  /// Run the physics list for a given material and particle.
  ///
  /// @param[in]     generator is the random number generator
  /// @param[in]     slab      is the passed material
  /// @param[in,out] particle  is the particle being updated
  /// @param[out]    generated is the container of generated particles
  /// @return Break condition, i.e. whether a process stoped the propagation
  ///
  /// @tparam generator_t must be a RandomNumberEngine
  template <typename generator_t>
  bool operator()(generator_t& generator, const Acts::MaterialProperties& slab,
                  Particle& particle, std::vector<Particle>& generated) const {
    static_assert(
        (true && ... &&
         std::is_same_v<bool, decltype(processes_t()(generator, slab, particle,
                                                     generated))>),
        "Not all processes conform to the expected interface");

    return impl(std::index_sequence_for<processes_t...>(), generator, slab,
                particle, generated);
  }

  /// Access a specific process by index.
  template <size_t I>
  std::tuple_element_t<I, std::tuple<processes_t...>>& get() {
    return std::get<I>(m_processes);
  }
  /// Access a specific process by type.
  template <typename process_t>
  process_t& get() {
    return std::get<process_t>(m_processes);
  }

  /// Disable a specific process by index.
  void disable(std::size_t i) { m_mask.set(i, true); }
  /// Disable a specific process by type.
  ///
  /// @warning Disables only the first of multiple processes of the same type.
  template <typename process_t>
  void disable() {
    return disable(Index<process_t, std::tuple<processes_t...>>::value);
  }

 private:
  // utility struct to retrieve index of the first matching type in the tuple.
  // from https://stackoverflow.com/a/18063608.
  template <class T, class Tuple>
  struct Index;
  template <class T, class... Types>
  struct Index<T, std::tuple<T, Types...>> {
    static constexpr std::size_t value = 0u;
  };
  template <class T, class U, class... Types>
  struct Index<T, std::tuple<U, Types...>> {
    static constexpr std::size_t value =
        1u + Index<T, std::tuple<Types...>>::value;
  };

  std::bitset<sizeof...(processes_t)> m_mask;
  std::tuple<processes_t...> m_processes;

  // compile-time index-based recursive function call for all processes
  template <typename generator_t>
  bool impl(std::index_sequence<>, generator_t&,
            const Acts::MaterialProperties&, Particle&,
            std::vector<Particle>&) const {
    return false;
  }
  template <std::size_t I0, std::size_t... INs, typename generator_t>
  bool impl(std::index_sequence<I0, INs...>, generator_t& generator,
            const Acts::MaterialProperties& slab, Particle& particle,
            std::vector<Particle>& generated) const {
    // only call process if it is not masked
    if (not m_mask[I0] and
        std::get<I0>(m_processes)(generator, slab, particle, generated)) {
      // exit early in case the process signals an abort
      return true;
    } else {
      return impl(std::index_sequence<INs...>(), generator, slab, particle,
                  generated);
    }
  }
};

}  // namespace ActsFatras
