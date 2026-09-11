// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges/ranges.hpp"

// Detray test include(s)
#include "detray/test/common/event_generator/random_track_generator_config.hpp"

// System include(s)
#include <algorithm>
#include <limits>
#include <memory>
#include <random>

namespace detray {

/// @brief Generates track states with random momentum directions.
///
/// Generates the phi and theta angles of the track momentum according to a
/// given random number distribution.
///
/// @tparam track_t the type of track parametrization that should be used.
/// @tparam generator_t source of random numbers
///
/// @note Since the random number generator might not be copy constructible,
/// neither is this generator. The iterators hold a reference to the rand
/// generator, which must not be invalidated during the iteration.
/// @note the random numbers are clamped to fit the phi/theta ranges. This can
/// effect distribution mean etc.
template <typename track_t,
          typename generator_t =
              detail::random_numbers<dscalar<typename track_t::algebra_type>>>
class random_track_generator
    : public detray::ranges::view_interface<
          random_track_generator<track_t, generator_t>> {
  using algebra_t = typename track_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

 public:
  using track_type = track_t;

  /// Configure how tracks are generated
  using configuration = random_track_generator_config<scalar_t>;

 private:
  /// @brief Nested iterator type that generates track states.
  struct iterator {
    using difference_type = std::ptrdiff_t;
    using value_type = track_t;
    using pointer = track_t*;
    using reference = track_t&;
    using iterator_category = detray::ranges::input_iterator_tag;

    constexpr iterator() = delete;

    DETRAY_HOST_DEVICE
    iterator(std::shared_ptr<generator_t> rand_gen, configuration cfg,
             std::size_t n_tracks_)
        : m_rnd_numbers{std::move(rand_gen)}, m_tracks{n_tracks_}, m_cfg{cfg} {}

    /// @returns whether we reached the end of iteration
    DETRAY_HOST_DEVICE
    constexpr bool operator==(const iterator& rhs) const {
      return rhs.m_tracks == m_tracks;
    }

    /// @returns the generator at its next position (prefix)
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> iterator& {
      ++m_tracks;
      return *this;
    }

    /// @returns the generator at its next position (postfix)
    DETRAY_HOST_DEVICE
    constexpr auto operator++(int) -> iterator& {
      auto tmp(m_tracks);
      ++m_tracks;
      return tmp;
    }

    /// @returns a track instance from random-generated momentum
    DETRAY_HOST_DEVICE
    track_t operator*() const {
      if (!m_rnd_numbers) {
        std::string err_str{"Invalid random number generator"};
        DETRAY_FATAL_HOST(err_str);
        throw std::invalid_argument(err_str);
      }

      const auto& ori = m_cfg.origin();
      const auto& ori_stddev = m_cfg.origin_stddev();

      const point3_t vtx =
          m_cfg.do_vertex_smearing()
              ? point3_t{m_rnd_numbers->normal(ori[0], ori_stddev[0]),
                         m_rnd_numbers->normal(ori[1], ori_stddev[1]),
                         m_rnd_numbers->normal(ori[2], ori_stddev[2])}
              : point3_t{ori[0], ori[1], ori[2]};

      scalar_t p_mag{(*m_rnd_numbers)(m_cfg.mom_range())};
      scalar_t phi{(*m_rnd_numbers)(m_cfg.phi_range())};
      scalar_t theta{(*m_rnd_numbers)(m_cfg.theta_range())};
      scalar_t sin_theta{math::sin(theta)};

      // Momentum direction from angles
      vector3_t mom{math::cos(phi) * sin_theta, math::sin(phi) * sin_theta,
                    math::cos(theta)};

      if constexpr (std::is_same_v<track_t, detail::ray<algebra_t>>) {
        mom = vector::normalize(mom);
      } else {
        sin_theta = (sin_theta == scalar_t{0.f})
                        ? std::numeric_limits<scalar_t>::epsilon()
                        : sin_theta;
        mom = (m_cfg.is_pT() ? 1.f / sin_theta : 1.f) * p_mag *
              vector::normalize(mom);
      }

      // Randomly flip the charge sign
      darray<double, 2> signs{1., -1.};
      const auto sign{static_cast<scalar_t>(
          signs[m_cfg.randomize_charge() ? m_rnd_numbers->coin_toss() : 0u])};

      return track_t{vtx, m_cfg.time(), mom, sign * m_cfg.charge()};
    }

    /// Random number generator
    std::shared_ptr<generator_t> m_rnd_numbers;

    /// How many tracks will be generated
    std::size_t m_tracks{0u};

    /// Configuration
    configuration m_cfg{};
  };

  std::shared_ptr<generator_t> m_gen{nullptr};
  configuration m_cfg{};

 public:
  using iterator_t = iterator;

  /// Default constructor
  constexpr random_track_generator() = default;

  /// Construct from external configuration
  DETRAY_HOST_DEVICE
  explicit constexpr random_track_generator(const configuration& cfg)
      : m_gen{std::make_shared<generator_t>(cfg.seed())}, m_cfg(cfg) {}

  /// Paramtetrized constructor for quick construction of simple tasks
  ///
  /// @note For more complex tasks, use the @c configuration type
  ///
  /// @param n_tracks_ the number of steps in the theta space
  /// @param mom_range the range of the track momentum (in GeV)
  /// @param charge charge of particle (e)
  DETRAY_HOST_DEVICE
  explicit random_track_generator(
      std::size_t n_tracks_,
      darray<scalar_t, 2> mom_range = {1.f * unit<scalar_t>::GeV,
                                       1.f * unit<scalar_t>::GeV},
      scalar_t charge = -1.f * unit<scalar_t>::e)
      : m_gen{std::make_shared<generator_t>()} {
    m_cfg.n_tracks(n_tracks_);
    m_cfg.mom_range(mom_range);
    m_cfg.charge(charge);
  }

  /// Move constructor
  DETRAY_HOST_DEVICE
  random_track_generator(random_track_generator&& other) noexcept
      : m_gen(std::move(other.m_gen)), m_cfg(std::move(other.m_cfg)) {}

  /// Access the configuration
  DETRAY_HOST_DEVICE
  constexpr configuration& config() { return m_cfg; }

  /// @returns the generator in initial state.
  /// @note the underlying random number generator has deleted copy
  /// constructor, so the iterator needs to be built from scratch
  DETRAY_HOST_DEVICE
  auto begin() noexcept -> iterator { return {m_gen, m_cfg, 0u}; }

  /// @returns the generator in end state
  DETRAY_HOST_DEVICE
  auto end() noexcept -> iterator { return {m_gen, m_cfg, m_cfg.n_tracks()}; }

  /// @returns the number of tracks that will be generated
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> std::size_t {
    return m_cfg.n_tracks();
  }
};

}  // namespace detray
