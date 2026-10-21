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
#include "detray/test/common/event_generator/uniform_track_generator_config.hpp"

// System include(s)
#include <algorithm>
#include <limits>
#include <memory>

namespace detray {

/// @brief Generates track states with momentum directions in a uniform angle
/// space.
///
/// It generates the track instances on the fly according to given parameters
/// and with the momentum direction determined by phi and theta angles, which
/// are advanced as the iteration proceeds. The angle space spans
/// theta [0, pi) x phi [-pi, pi), while the step sizes (and with them
/// the number of generated tracks) are configurable.
///
/// @tparam track_t the type of track parametrization that should be used.
template <typename track_t,
          typename generator_t =
              detail::random_numbers<dscalar<typename track_t::algebra_type>>>
class uniform_track_generator
    : public detray::ranges::view_interface<
          uniform_track_generator<track_t, generator_t>> {
  using algebra_t = typename track_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

 public:
  using track_type = track_t;

  /// Configure how tracks are generated
  using configuration = uniform_track_generator_config<scalar_t>;

 private:
  /// @brief Nested iterator type that generates track states.
  struct iterator {
    using difference_type = std::ptrdiff_t;
    using value_type = track_t;
    using pointer = track_t*;
    using reference = track_t&;
    using iterator_category = detray::ranges::input_iterator_tag;

    constexpr iterator() = default;

    DETRAY_HOST_DEVICE
    constexpr iterator(std::shared_ptr<generator_t> rand_gen, configuration cfg,
                       std::size_t iph = 1u, std::size_t ith = 0u)
        : m_rnd_numbers{std::move(rand_gen)},
          m_cfg{cfg},
          m_phi_step_size{(cfg.phi_range()[1] - cfg.phi_range()[0]) /
                          static_cast<scalar_t>(cfg.phi_steps())},
          m_theta_step_size{(cfg.theta_range()[1] - cfg.theta_range()[0]) /
                            static_cast<scalar_t>(cfg.theta_steps() - 1u)},
          m_eta_step_size{(cfg.eta_range()[1] - cfg.eta_range()[0]) /
                          static_cast<scalar_t>(cfg.eta_steps() - 1u)},
          m_phi{cfg.phi_range()[0]},
          m_theta{cfg.uniform_eta() ? get_theta(cfg.eta_range()[0])
                                    : cfg.theta_range()[0]},
          i_phi{iph},
          i_theta{ith} {}

    /// @returns whether we reached end of angle space
    DETRAY_HOST_DEVICE
    constexpr bool operator==(const iterator& rhs) const {
      return rhs.i_phi == i_phi && rhs.i_theta == i_theta;
    }

    /// Iterate through angle space according to given step sizes.
    ///
    /// @returns the generator at its next position.
    DETRAY_HOST_DEVICE
    constexpr auto operator++() -> iterator& {
      // Check theta range according to step size
      if (i_theta < m_cfg.theta_steps()) {
        // Check phi sub-range
        if (i_phi < m_cfg.phi_steps()) {
          // Calculate new phi in the given range
          m_phi = m_cfg.phi_range()[0] +
                  static_cast<scalar_t>(i_phi) * m_phi_step_size;
          ++i_phi;
          return *this;
        }
        // Reset phi range
        i_phi = 1;
        m_phi = m_cfg.phi_range()[0];

        // Calculate new theta in the given range
        ++i_theta;

        if (m_cfg.uniform_eta()) {
          const scalar_t eta = m_cfg.eta_range()[0] +
                               static_cast<scalar_t>(i_theta) * m_eta_step_size;
          m_theta = get_theta(eta);
        } else {
          m_theta = m_cfg.theta_range()[0] +
                    static_cast<scalar_t>(i_theta) * m_theta_step_size;
        }
      }
      return *this;
    }

    /// @returns the generator at its next position (postfix)
    DETRAY_HOST_DEVICE
    constexpr auto operator++(int) -> iterator& {
      auto tmp(*this);
      ++(*this);
      return tmp;
    }

    /// @returns a track instance from generated momentum direction
    DETRAY_HOST_DEVICE
    track_t operator*() const {
      if (!m_rnd_numbers) {
        std::string err_str{"Invalid random number generator"};
        DETRAY_FATAL_HOST(err_str);
        throw std::invalid_argument(err_str);
      }

      scalar_t sin_theta{math::sin(m_theta)};

      // Momentum direction from angles
      vector3_t p{math::cos(m_phi) * sin_theta, math::sin(m_phi) * sin_theta,
                  math::cos(m_theta)};

      // Magnitude of momentum
      if constexpr (std::is_same_v<track_t, detail::ray<algebra_t>>) {
        p = vector::normalize(p);
      } else {
        sin_theta = (sin_theta == scalar_t{0.f})
                        ? std::numeric_limits<scalar_t>::epsilon()
                        : sin_theta;
        p = (m_cfg.is_pT() ? 1.f / sin_theta : 1.f) * m_cfg.m_p_mag *
            vector::normalize(p);
      }

      const auto& ori = m_cfg.origin();

      // Randomly flip the charge sign
      darray<double, 2> signs{1., -1.};
      const auto sign{static_cast<scalar_t>(
          signs[m_cfg.randomize_charge() ? m_rnd_numbers->coin_toss() : 0u])};

      return track_t{
          {ori[0], ori[1], ori[2]}, m_cfg.time(), p, sign * m_cfg.charge()};
    }

    /// Random number generator
    std::shared_ptr<generator_t> m_rnd_numbers;

    /// Current configuration
    configuration m_cfg{};

    /// Angular step sizes
    scalar_t m_phi_step_size{0.f};
    scalar_t m_theta_step_size{0.f};
    scalar_t m_eta_step_size{0.f};

    /// Phi and theta angles of momentum direction
    scalar_t m_phi{-constant<scalar_t>::pi};
    scalar_t m_theta{0.f};

    /// Iteration indices
    std::size_t i_phi{0u};
    std::size_t i_theta{0u};

   private:
    /// @returns the theta angle for a given @param eta value
    DETRAY_HOST_DEVICE
    scalar_t get_theta(const scalar_t eta) {
      return 2.f * math::atan(math::exp(-eta));
    }
  };

  std::shared_ptr<generator_t> m_gen{
      std::make_shared<generator_t>(configuration{}.seed())};
  configuration m_cfg{};

 public:
  using iterator_t = iterator;

  /// Default constructor
  constexpr uniform_track_generator() = default;

  /// Construct from external configuration @param cfg
  DETRAY_HOST_DEVICE
  explicit constexpr uniform_track_generator(configuration cfg)
      : m_gen{std::make_shared<generator_t>(cfg.seed())}, m_cfg{cfg} {}

  /// Paramtetrized constructor for quick construction of simple tasks
  ///
  /// @note For more complex tasks, use the @c configuration type
  ///
  /// @param n_theta the number of steps in the theta space
  /// @param n_phi the number of steps in the phi space
  /// @param p_mag magnitude of the track momentum (in GeV)
  /// @param uniform_eta uniformly step through eta space instead of theta
  /// @param charge charge of particle (e)
  DETRAY_HOST_DEVICE
  uniform_track_generator(std::size_t n_phi, std::size_t n_theta,
                          scalar_t p_mag = 1.f * unit<scalar_t>::GeV,
                          bool uniform_eta = false,
                          scalar_t charge = -1.f * unit<scalar_t>::e)
      : m_gen{std::make_shared<generator_t>()} {
    m_cfg.phi_steps(n_phi).theta_steps(n_theta);
    m_cfg.uniform_eta(uniform_eta);
    m_cfg.p_tot(p_mag);
    m_cfg.charge(charge);
  }

  /// Move constructor
  DETRAY_HOST_DEVICE
  uniform_track_generator(uniform_track_generator&& other) noexcept
      : m_gen(std::move(other.m_gen)), m_cfg(std::move(other.m_cfg)) {}

  /// Copy assignment operator
  DETRAY_HOST_DEVICE
  constexpr uniform_track_generator& operator=(
      const uniform_track_generator& other) {
    m_cfg = other.m_cfg;
    m_gen = other.m_gen;
    return *this;
  }

  /// Access the configuration
  DETRAY_HOST_DEVICE
  constexpr configuration& config() { return m_cfg; }

  /// @returns the generator in initial state: Default values reflect the
  /// first phi angle iteration.
  DETRAY_HOST_DEVICE
  constexpr auto begin() noexcept -> iterator { return {m_gen, m_cfg, 1u, 0u}; }

  /// @returns the generator in end state
  DETRAY_HOST_DEVICE
  constexpr auto end() noexcept -> iterator {
    return {m_gen, m_cfg, 1u, m_cfg.theta_steps()};
  }

  /// @returns the number of tracks that will be generated
  DETRAY_HOST_DEVICE
  constexpr auto size() const noexcept -> std::size_t {
    return m_cfg.phi_steps() * m_cfg.theta_steps();
  }
};

}  // namespace detray
