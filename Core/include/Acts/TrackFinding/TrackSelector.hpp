// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <ostream>
#include <vector>

namespace Acts {

/// Class which performs filtering of tracks. It accepts an input and an output
/// track container and uses the built-in copy facility to copy tracks into the
/// output container.
class TrackSelector {
  static constexpr double inf = std::numeric_limits<double>::infinity();

 public:
  /// Configuration of a set of cuts for a single eta bin
  /// Default construction yields a set of cuts that accepts everything.
  struct CutConfig {
    // Minimum/maximum local positions.
    double loc0Min = -inf;
    double loc0Max = inf;
    double loc1Min = -inf;
    double loc1Max = inf;
    // Minimum/maximum track time.
    double timeMin = -inf;
    double timeMax = inf;
    // Direction cuts.
    double phiMin = -inf;
    double phiMax = inf;
    double etaMin = -inf;
    double etaMax = inf;
    double absEtaMin = 0.0;
    double absEtaMax = inf;
    // Momentum cuts.
    double ptMin = 0.0;
    double ptMax = inf;

    std::size_t minMeasurements = 0;

    // Helper factory functions to produce a populated config object more
    // conveniently

    /// Set loc0 acceptance range
    /// @param min Minimum value
    /// @param max Maximum value
    /// @return Reference to this object
    CutConfig& loc0(double min, double max);

    /// Set loc1 acceptance range
    /// @param min Minimum value
    /// @param max Maximum value
    /// @return Reference to this object
    CutConfig& loc1(double min, double max);

    /// Set time acceptance range
    /// @param min Minimum value
    /// @param max Maximum value
    /// @return Reference to this object
    CutConfig& time(double min, double max);

    /// Set phi acceptance range
    /// @param min Minimum value
    /// @param max Maximum value
    /// @return Reference to this object
    CutConfig& phi(double min, double max);

    /// Set the eta acceptance range
    /// @param min Minimum value
    /// @param max Maximum value
    /// @return Reference to this object
    CutConfig& eta(double min, double max);

    /// Set the absolute eta acceptance range
    /// @param min Minimum value
    /// @param max Maximum value
    /// @return Reference to this object
    CutConfig& absEta(double min, double max);

    /// Set the pt acceptance range
    /// @param min Minimum value
    /// @param max Maximum value
    /// @return Reference to this object
    CutConfig& pt(double min, double max);

    /// Print this set of cuts to an output stream
    /// @param os Output stream
    /// @param cuts Cuts to print
    /// @return Reference to the output stream
    friend std::ostream& operator<<(std::ostream& os, const CutConfig& cuts);
  };

  /// Main config object for the track selector. Combines a set of cut
  /// configurations and corresponding eta bins
  struct Config {
    /// Cut sets for each eta bin
    std::vector<CutConfig> cutSets = {};

    /// Eta bin edges for varying cuts by eta
    std::vector<double> absEtaEdges = {};

    /// Get the number of eta bins
    /// @return Number of eta bins
    std::size_t nEtaBins() const { return absEtaEdges.size() - 1; }

    /// Construct an empty (accepts everything) configuration.
    /// Results in a single cut set and one abs eta bin from 0 to infinity.
    Config() : cutSets{{}}, absEtaEdges{{0, inf}} {};

    /// Constructor to create a config object that is not upper-bounded.
    /// This is useful to use the "fluent" API to populate the configuration.
    /// @param etaMin Minimum eta bin edge
    Config(double etaMin) : cutSets{}, absEtaEdges{etaMin} {}

    /// Constructor from a vector of eta bin edges. This automatically
    /// initializes all the cuts to be the same for all eta and be essentially
    /// no-op.
    /// @param absEtaEdgesIn is the vector of eta bin edges
    Config(std::vector<double> absEtaEdgesIn)
        : absEtaEdges{std::move(absEtaEdgesIn)} {
      cutSets.resize(absEtaEdges.size() - 1);
    }

    /// Auto-converting constructor from a single cut configuration.
    /// Results in a single absolute eta bin from 0 to infinity.
    Config(CutConfig cutSet) : cutSets{cutSet}, absEtaEdges{{0, inf}} {}

    /// Add a new eta bin with the given upper bound.
    /// @param etaMax Upper bound of the new eta bin
    /// @param callback Callback to configure the cuts for this eta bin
    /// @return Reference to this object
    Config& addCuts(double etaMax,
                    const std::function<void(CutConfig&)>& callback = {});

    /// Add a new eta bin with an upper bound of +infinity.
    /// @param callback Callback to configure the cuts for this eta bin
    /// @return Reference to this object
    Config& addCuts(const std::function<void(CutConfig&)>& callback = {});

    /// Print this configuration to an output stream
    /// @param os Output stream
    /// @param cfg Configuration to print
    /// @return Reference to the output stream
    friend std::ostream& operator<<(std::ostream& os, const Config& cfg);

    /// Get the index of the eta bin for a given eta
    /// @param eta Eta value
    /// @return Index of the eta bin
    std::size_t binIndex(double eta) const;

    /// Get the cuts for a given eta
    /// @param eta Eta value
    /// @return Cuts for the given eta
    const CutConfig& getCuts(double eta) const;
  };

  /// Constructor from a config object
  /// @param config is the configuration object
  TrackSelector(const Config& config);

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

inline TrackSelector::CutConfig& TrackSelector::CutConfig::loc0(double min,
                                                                double max) {
  loc0Min = min;
  loc0Max = max;
  return *this;
}

inline TrackSelector::CutConfig& TrackSelector::CutConfig::loc1(double min,
                                                                double max) {
  loc1Min = min;
  loc1Max = max;
  return *this;
}

inline TrackSelector::CutConfig& TrackSelector::CutConfig::time(double min,
                                                                double max) {
  timeMin = min;
  timeMax = max;
  return *this;
}

inline TrackSelector::CutConfig& TrackSelector::CutConfig::phi(double min,
                                                               double max) {
  phiMin = min;
  phiMax = max;
  return *this;
}

inline TrackSelector::CutConfig& TrackSelector::CutConfig::eta(double min,
                                                               double max) {
  if (absEtaMin != 0.0 || absEtaMax != inf) {
    throw std::invalid_argument(
        "Cannot set both eta and absEta cuts in the same cut set");
  }
  etaMin = min;
  etaMax = max;
  return *this;
}

inline TrackSelector::CutConfig& TrackSelector::CutConfig::absEta(double min,
                                                                  double max) {
  if (etaMin != -inf || etaMax != inf) {
    throw std::invalid_argument(
        "Cannot set both eta and absEta cuts in the same cut set");
  }
  absEtaMin = min;
  absEtaMax = max;
  return *this;
}

inline TrackSelector::CutConfig& TrackSelector::CutConfig::pt(double min,
                                                              double max) {
  ptMin = min;
  ptMax = max;
  return *this;
}

inline std::ostream& operator<<(std::ostream& os,
                                const TrackSelector::CutConfig& cuts) {
  auto print = [&](const char* name, const auto& min, const auto& max) {
    os << " - " << min << " <= " << name << " < " << max << "\n";
  };

  print("loc0", cuts.loc0Min, cuts.loc0Max);
  print("loc1", cuts.loc1Min, cuts.loc1Max);
  print("time", cuts.timeMin, cuts.timeMax);
  print("phi", cuts.phiMin, cuts.phiMax);
  print("eta", cuts.etaMin, cuts.etaMax);
  print("absEta", cuts.absEtaMin, cuts.absEtaMax);
  print("pt", cuts.ptMin, cuts.ptMax);
  os << " - " << cuts.minMeasurements << " <= nMeasurements\n";

  return os;
}

inline TrackSelector::Config& TrackSelector::Config::addCuts(
    double etaMax, const std::function<void(CutConfig&)>& callback) {
  if (etaMax <= absEtaEdges.back()) {
    throw std::invalid_argument{
        "Abs Eta bin edges must be in increasing order"};
  }

  if (etaMax < 0.0) {
    throw std::invalid_argument{"Abs Eta bin edges must be positive"};
  }

  absEtaEdges.push_back(etaMax);
  cutSets.emplace_back();
  if (callback) {
    callback(cutSets.back());
  }
  return *this;
}

inline TrackSelector::Config& TrackSelector::Config::addCuts(
    const std::function<void(CutConfig&)>& callback) {
  return addCuts(inf, callback);
}

inline std::size_t TrackSelector::Config::binIndex(double eta) const {
  if (std::abs(eta) >= absEtaEdges.back()) {
    throw std::invalid_argument{"Eta is outside the abs eta bin edges"};
  }

  auto binIt =
      std::upper_bound(absEtaEdges.begin(), absEtaEdges.end(), std::abs(eta));
  std::size_t index = std::distance(absEtaEdges.begin(), binIt) - 1;
  return index;
}

inline const TrackSelector::CutConfig& TrackSelector::Config::getCuts(
    double eta) const {
  return nEtaBins() == 1 ? cutSets[0] : cutSets[binIndex(eta)];
}

inline std::ostream& operator<<(std::ostream& os,
                                const TrackSelector::Config& cfg) {
  os << "TrackSelector::Config:\n";

  for (std::size_t i = 1; i < cfg.absEtaEdges.size(); i++) {
    os << cfg.absEtaEdges[i - 1] << " <= eta < " << cfg.absEtaEdges[i] << "\n";
    os << cfg.cutSets[i - 1];
  }

  return os;
}

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
  const auto absEta = std::abs(eta);

  if (absEta < m_cfg.absEtaEdges.front() ||
      absEta >= m_cfg.absEtaEdges.back()) {
    return false;
  }

  const CutConfig& cuts = m_cfg.getCuts(eta);

  return within(track.transverseMomentum(), cuts.ptMin, cuts.ptMax) and
         within(std::abs(eta), cuts.absEtaMin, cuts.absEtaMax) and
         within(eta, cuts.etaMin, cuts.etaMax) and
         within(track.phi(), cuts.phiMin, cuts.phiMax) and
         within(track.loc0(), cuts.loc0Min, cuts.loc0Max) and
         within(track.loc1(), cuts.loc1Min, cuts.loc1Max) and
         within(track.time(), cuts.timeMin, cuts.timeMax) and
         checkMin(track.nMeasurements(), cuts.minMeasurements);
}

inline TrackSelector::TrackSelector(const TrackSelector::Config& config)
    : m_cfg(config) {
  if (m_cfg.cutSets.size() != m_cfg.nEtaBins()) {
    throw std::invalid_argument{
        "TrackSelector cut / eta bin configuration is inconsistent"};
  }

  if (m_cfg.nEtaBins() == 1) {
    static const std::vector<double> infVec = {0, inf};
    bool limitEta = m_cfg.absEtaEdges != infVec;

    const CutConfig& cuts = m_cfg.cutSets[0];

    if (limitEta && (cuts.etaMin != -inf || cuts.etaMax != inf ||
                     cuts.absEtaMin != 0.0 || cuts.absEtaMax != inf)) {
      throw std::invalid_argument{
          "Explicit eta cuts are only valid for single eta bin"};
    }
  }
}

}  // namespace Acts
