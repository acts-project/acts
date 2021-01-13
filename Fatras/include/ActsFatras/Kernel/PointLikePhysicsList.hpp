// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/EventData/Particle.hpp"

#include <bitset>
#include <cassert>
#include <limits>
#include <tuple>

namespace ActsFatras {

/// List of point-like physics processes.
template <typename... processes_t>
class PointLikePhysicsList {
 public:
  struct Selection {
    Particle::Scalar x0Limit =
        std::numeric_limits<Particle::Scalar>::infinity();
    Particle::Scalar l0Limit =
        std::numeric_limits<Particle::Scalar>::infinity();
    size_t x0Process = SIZE_MAX;
    size_t l0Process = SIZE_MAX;
  };

  /// Disable a specific process by type.
  template <typename process_t>
  void disable() {
    m_mask.set(Index<process_t, Processes>::value);
  }

  /// Access a specific process by type.
  template <typename process_t>
  process_t& get() {
    return std::get<process_t>(m_processes);
  }

  /// Arm the physics lists by generating limits and selecting processes.
  ///
  /// @tparam generator_t must be a RandomNumberEngine
  /// @param[in] rng is the random number generator
  /// @param[in] particle is the initial particle state
  /// @return X0/L0 limits for the particle and the process index that should be
  ///   executed once the limit has been reached.
  template <typename generator_t>
  Selection arm(generator_t& rng, const Particle& particle) const {
    Selection selection;
    armImpl(rng, particle, selection,
            std::index_sequence_for<processes_t...>());
    return selection;
  }

  /// Run one physics processes for the particle.
  ///
  /// @tparam generator_t must be a RandomNumberEngine
  /// @param[in] rng is the random number generator
  /// @param[in] process is the index of the process to be executed
  /// @param[in,out] particle  is the particle being updated
  /// @param[out] generated is the container of generated particles
  /// @return Break condition, i.e. whether a process killed the particle
  ///
  /// The process index is expected to originate from a previous `arm(...)`
  /// call, but this is not enforced. How to select the correct process requires
  /// more information that is not available here.
  template <typename generator_t>
  bool run(generator_t& rng, size_t process, Particle& particle,
           std::vector<Particle>& generated) const {
    static_assert((true && ... &&
                   std::is_same_v<bool, decltype(processes_t().run(
                                            rng, particle, generated))>),
                  "Not all processes conform to the expected interface");

    return runImpl(rng, process, particle, generated,
                   std::index_sequence_for<processes_t...>());
  }

 private:
  // TODO check that all processes are unique types.

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

  using Mask = std::bitset<sizeof...(processes_t)>;
  using Processes = std::tuple<processes_t...>;

  // allow processes to be masked. defaults to zeros -> no masked processes
  Mask m_mask;
  Processes m_processes;

  // for the `arm` call, we need to iterate over all available processes and
  // select the ones that generate the smallest limits. this is done using an
  // index-based compile-time recursive call.
  template <typename generator_t>
  void armImpl(generator_t&, const Particle&, Selection&,
               std::index_sequence<>) const {}
  template <typename generator_t, std::size_t I0, std::size_t... INs>
  void armImpl(generator_t& rng, const Particle& particle, Selection& selection,
               std::index_sequence<I0, INs...>) const {
    // only arm the process if it is not masked
    if (not m_mask[I0]) {
      auto [x0Limit, l0Limit] =
          std::get<I0>(m_processes).generatePathLimits(rng, particle);
      if (x0Limit < selection.x0Limit) {
        selection.x0Limit = x0Limit;
        selection.x0Process = I0;
      }
      if (l0Limit < selection.l0Limit) {
        selection.l0Limit = l0Limit;
        selection.l0Process = I0;
      }
    }
    // continue with the remaining processes
    armImpl(rng, particle, selection, std::index_sequence<INs...>());
  }

  // for the `run` call we need call just one process. since we can not select a
  // tuple element with a runtime index, we need to iterate over all processes
  // with a compile-time recursive function until we reach the requested one.
  template <typename generator_t>
  bool runImpl(generator_t&, size_t, Particle&, std::vector<Particle>&,
               std::index_sequence<>) const {
    // the requested process index is outside the possible range. **do not**
    // treat this as an error to simplify the case of an empty physics lists or
    // a default process index.
    return false;
  }
  template <typename generator_t, size_t I0, size_t... INs>
  bool runImpl(generator_t& rng, size_t process, Particle& particle,
               std::vector<Particle>& generated,
               std::index_sequence<I0, INs...>) const {
    if (I0 == process) {
      if (m_mask[I0]) {
        // the selected process is masked. since nothing is executed the
        // particle continues to be alive; not a break condition.
        return false;
      }
      return std::get<I0>(m_processes).run(rng, particle, generated);
    }
    // continue the iteration with the remaining processes
    return runImpl(rng, process, particle, generated,
                   std::index_sequence<INs...>());
  }
};

}  // namespace ActsFatras
