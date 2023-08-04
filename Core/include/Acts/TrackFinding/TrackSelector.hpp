// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>

namespace Acts {

/// Class which performs filtering of tracks. It accepts an input and an output
/// track container and uses the built-in copy facility to copy tracks into the
/// output container.
class TrackSelector {
 public:
  struct Config {
    // Minimum/maximum local positions.
    double loc0Min = -std::numeric_limits<double>::infinity();
    double loc0Max = std::numeric_limits<double>::infinity();
    double loc1Min = -std::numeric_limits<double>::infinity();
    double loc1Max = std::numeric_limits<double>::infinity();
    // Minimum/maximum track time.
    double timeMin = -std::numeric_limits<double>::infinity();
    double timeMax = std::numeric_limits<double>::infinity();
    // Direction cuts.
    double phiMin = -std::numeric_limits<double>::infinity();
    double phiMax = std::numeric_limits<double>::infinity();
    double etaMin = -std::numeric_limits<double>::infinity();
    double etaMax = std::numeric_limits<double>::infinity();
    double absEtaMin = 0.0;
    double absEtaMax = std::numeric_limits<double>::infinity();
    // Momentum cuts.
    double ptMin = 0.0;
    double ptMax = std::numeric_limits<double>::infinity();

    std::size_t minMeasurements = 0;
  };

  /// Constructor from a config object
  /// @param config is the configuration object
  TrackSelector(const Config& config) : m_cfg(config) {}

  /// Select tracks from an input container and copy them into an output
  /// container
  /// @tparam input_tracks_t is the type of the input track container
  /// @tparam output_tracks_t is the type of the output track container
  /// @param inputTracks is the input track container
  /// @param outputTracks is the output track container
  template <typename input_tracks_t, typename output_tracks_t>
  void selectTracks(const input_tracks_t& inputTracks,
                    output_tracks_t& outputTracks) const;

  /// Helper function to check if a track is valid
  /// @tparam track_proxy_t is the type of the track proxy
  /// @param track is the track proxy
  /// @return true if the track is valid
  template <typename track_proxy_t>
  bool isValidTrack(const track_proxy_t& track) const;

  /// Get readonly access to the config parameters
  /// @return the config object
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

template <typename input_tracks_t, typename output_tracks_t>
void TrackSelector::selectTracks(const input_tracks_t& inputTracks,
                                 output_tracks_t& outputTracks) const {
  for (auto track : inputTracks) {
    if (!isValidTrack(track)) {
      continue;
    }
    auto destProxy = outputTracks.getTrack(outputTracks.addTrack());
    destProxy.copyFrom(track, false);
    destProxy.tipIndex() = track.tipIndex();
  }
}

template <typename track_proxy_t>
bool TrackSelector::isValidTrack(const track_proxy_t& track) const {
  auto checkMin = [](auto x, auto min) { return min <= x; };
  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  const auto theta = track.theta();
  const auto eta = -std::log(std::tan(theta / 2));
  return within(track.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) and
         within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) and
         within(eta, m_cfg.etaMin, m_cfg.etaMax) and
         within(track.phi(), m_cfg.phiMin, m_cfg.phiMax) and
         within(track.loc0(), m_cfg.loc0Min, m_cfg.loc0Max) and
         within(track.loc1(), m_cfg.loc1Min, m_cfg.loc1Max) and
         within(track.time(), m_cfg.timeMin, m_cfg.timeMax) and
         checkMin(track.nMeasurements(), m_cfg.minMeasurements);
}

}  // namespace Acts
