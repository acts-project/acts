/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// traccc include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/utils/particle.hpp"

// detray include(s).
#include <detray/propagator/propagation_config.hpp>

// System include(s)
#include <iostream>

namespace traccc {

/// Indicate the type of track state candidate collection
enum class smoother_type : std::uint_least8_t {
  e_none,
  e_mbf,
  e_kalman,
  e_unknown,
};

TRACCC_HOST inline std::ostream &operator<<(std::ostream &os,
                                            smoother_type st) {
  using enum smoother_type;
  switch (st) {
    case e_none:
      os << "none";
      break;
    case e_mbf:
      os << "MBF";
      break;
    case e_kalman:
      os << "Kalman";
      break;
    default:
      os << "unknown";
  }
  return os;
}

TRACCC_HOST inline std::istream &operator>>(std::istream &is,
                                            smoother_type &st) {
  using enum smoother_type;
  std::string tmp_str;
  is >> tmp_str;
  if (tmp_str == "none") {
    st = e_none;
  } else if (tmp_str == "MBF") {
    st = e_mbf;
  } else if (tmp_str == "Kalman") {
    st = e_kalman;
  } else {
    st = e_unknown;
  }
  return is;
}

/// Configuration struct for track finding
struct finding_config {
  /// Maxmimum number of branches per seed
  unsigned int max_num_branches_per_seed = 10;

  /// Maximum number of branches per surface
  unsigned int max_num_branches_per_surface = 1;

  /// Min/Max number of track candidates per track
  unsigned int min_track_candidates_per_track = 3;
  unsigned int max_track_candidates_per_track = 100;

  /// Maximum allowed number of skipped steps per candidate
  unsigned int max_num_skipping_per_cand = 3;

  /// Maximum number of consecutive holes
  unsigned int max_num_consecutive_skipped = 1;

  /// Maximum number of tracks per measurement; if zero, don't prune
  unsigned int max_num_tracks_per_measurement = 0;

  /// If `max_num_tracks_per_measurement` is enabled, i.e. if it is non-zero
  /// then this value determines the minimum fraction of measurements in
  /// each track that vote for it
  float min_measurement_voting_fraction = 0.5f;

  /// Run the progressive Kalman filter (PKF) for track finding if branching
  /// is turned off
  bool run_pkf = false;
  /// The type of smoother to be run in track finding
  smoother_type run_smoother = smoother_type::e_mbf;

  /// Minimum step length that track should make to reach the next surface. It
  /// should be set higher than the overstep tolerance not to make it stay on
  /// the same surface
  float min_step_length_for_next_surface = 1.2f * traccc::unit<float>::mm;
  /// Maximum step counts that track can make to reach the next surface
  unsigned int max_step_counts_for_next_surface = 100;

  /// Maximum Chi-square that is allowed for branching
  float chi2_max = 30.f;

  /// Propagation configuration
  detray::propagation::config propagation{};
  /// Measurement calibration configuration
  measurement_selector::config meas_calibration{};
  // Fitting config for the Kalman smoother
  traccc::fitting_config kalman_smoother;

  /// Minimum momentum for reconstructed tracks
  float min_p = 100.f * traccc::unit<float>::MeV;
  float min_pT = 600.f * traccc::unit<float>::MeV;

  /// Particle hypothesis
  traccc::pdg_particle<traccc::scalar> ptc_hypothesis =
      traccc::pion_plus<traccc::scalar>();

  /// Minimum track length in order for a track to be a candidate for
  /// duplicate removal.
  ///
  /// @warning This parameter should be greater than or equal to 3 under all
  /// circumstances!
  unsigned int duplicate_removal_minimum_length = 5u;

  /// @name Performance parameters
  /// These parameters impact only compute performance; any case in which a
  /// change in these parameters effects a change in _physics_ performance
  /// should be considered a bug.
  /// @{
  /// @brief The number of links to reserve space for, per seed.
  ///
  /// This parameter describes the number of links which we reserve per seed.
  /// If this number turns out to be too low, the track finding algorithm
  /// will automatically resize it, but this comes at the cost of
  /// performance.
  ///
  /// @note This parameter affects GPU-based track finding only.
  unsigned int initial_links_per_seed = 100;
  /// @}
};

}  // namespace traccc
