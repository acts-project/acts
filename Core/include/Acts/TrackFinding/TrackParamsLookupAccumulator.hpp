// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <concepts>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace {

/// @brief Shorthand for track parameters
template <class parameters_t>
concept TrackParameters = Acts::FreeTrackParametersConcept<parameters_t> ||
                          Acts::BoundTrackParametersConcept<parameters_t>;

/// @brief Shorthand for GenericBoundTrackParameters
template <class parameters_t>
concept IsGenericBound =
    std::same_as<parameters_t, Acts::GenericBoundTrackParameters<
                                   typename parameters_t::ParticleHypothesis>>;

template <typename T>
using remove_shared = std::remove_reference_t<decltype(*std::declval<T>())>;

/// @brief Concept that restricts the type of the
/// accumulation grid cell
template <typename grid_t>
concept TrackParamsGrid = requires {
  typename grid_t::value_type::first_type;
  typename grid_t::value_type::second_type;

  requires TrackParameters<
      remove_shared<typename grid_t::value_type::first_type>>;
  requires TrackParameters<
      remove_shared<typename grid_t::value_type::first_type>>;

  requires requires(typename grid_t::value_type val) {
    {
      val.first
    } -> std::same_as<std::shared_ptr<remove_shared<decltype(val.first)>>&>;
    {
      val.second
    } -> std::same_as<std::shared_ptr<remove_shared<decltype(val.first)>>&>;
    { val.second } -> std::same_as<decltype(val.first)&>;
  };
};

}  // namespace

namespace Acts {

/// @brief Class to accumulate and average track lookup tables
///
/// @tparam Grid type for track parameters accumulation
///
/// This class is used to accumulate track parameters in
/// reference layer grids and average them to create a lookup
/// table for track parameter estimation in seeding
///
/// @note Geometry context is left to be handled by the user
/// outside of accumulation
template <TrackParamsGrid grid_t>
class TrackParamsLookupAccumulator {
 public:
  using LookupGrid = grid_t;
  using TrackParameters = typename std::pointer_traits<
      typename grid_t::value_type::first_type>::element_type;

  /// @brief Constructor
  explicit TrackParamsLookupAccumulator(grid_t grid)
      : m_grid(std::move(grid)) {}

  /// @brief Add track parameters to the accumulator
  ///
  /// @param ipTrackParameters Track parameters at the IP
  /// @param refTrackParameters Track parameters at the reference layer
  /// @param position Local position of the track hit on the reference layer
  void addTrack(const TrackParameters& ipTrackParameters,
                const TrackParameters& refTrackParameters,
                const Vector2& position) {
    std::lock_guard<std::mutex> lock(m_gridMutex);

    auto bin = m_grid.localBinsFromPosition(position);

    if (m_countGrid[bin] == 0) {
      m_grid.atLocalBins(bin).first =
          std::make_shared<TrackParameters>(ipTrackParameters);
      m_grid.atLocalBins(bin).second =
          std::make_shared<TrackParameters>(refTrackParameters);

      m_countGrid.at(bin)++;
      return;
    }

    *m_grid.atLocalBins(bin).first =
        addTrackParameters(*m_grid.atLocalBins(bin).first, ipTrackParameters);
    *m_grid.atLocalBins(bin).second =
        addTrackParameters(*m_grid.atLocalBins(bin).second, refTrackParameters);
    m_countGrid.at(bin)++;
  }

  /// @brief Finalize the lookup table
  ///
  /// @return Grid with the bin track parameters averaged
  LookupGrid finalizeLookup() {
    auto meanTrack = [&](const TrackParameters& track, std::size_t count) {
      if constexpr (IsGenericBound<TrackParameters>) {
        Acts::GeometryContext gctx;

        auto res = TrackParameters::create(
            track.referenceSurface().getSharedPtr(), gctx,
            track.fourPosition(gctx) / count, track.momentum().normalized(),
            count * track.charge() / track.momentum().norm(),
            track.covariance(), track.particleHypothesis());

        if (!res.ok()) {
          throw std::invalid_argument("Bound track grid finalization failed");
        }
        return res.value();
      } else {
        return TrackParameters(track.fourPosition() / count,
                               track.momentum().normalized(),
                               count * track.charge() / track.momentum().norm(),
                               track.covariance(), track.particleHypothesis());
      }
    };

    for (auto [bin, count] : m_countGrid) {
      if (count == 0) {
        continue;
      }
      *m_grid.atLocalBins(bin).first =
          meanTrack(*m_grid.atLocalBins(bin).first, count);
      *m_grid.atLocalBins(bin).second =
          meanTrack(*m_grid.atLocalBins(bin).second, count);
    }

    return m_grid;
  }

 private:
  /// @brief Add two track parameters
  ///
  /// @param a First track parameter in the sum
  /// @param b Second track parameter in the sum
  ///
  /// @return Sum of track parameters a + b
  ///
  /// @note Covariances of the track parameters
  /// are not added and instead assumed to be
  /// generated by the same random process for
  /// both a and b, making its averaging redundant
  TrackParameters addTrackParameters(const TrackParameters& a,
                                     const TrackParameters& b) {
    if (a.particleHypothesis() != b.particleHypothesis()) {
      throw std::invalid_argument(
          "Cannot accumulate track parameters with different particle "
          "hypotheses");
    }
    if (a.charge() != b.charge()) {
      throw std::invalid_argument(
          "Cannot accumulate track parameters with different charges");
    }
    if constexpr (IsGenericBound<TrackParameters>) {
      if (a.referenceSurface() != b.referenceSurface()) {
        throw std::invalid_argument(
            "Cannot accumulate bound track parameters with different reference "
            "surfaces");
      }
    }

    Acts::Vector3 momentum = a.momentum() + b.momentum();

    // Assume track parameters being i.i.d.
    if constexpr (IsGenericBound<TrackParameters>) {
      Acts::GeometryContext gctx;

      Acts::Vector4 fourPosition = a.fourPosition(gctx) + b.fourPosition(gctx);

      auto res = TrackParameters::create(
          a.referenceSurface().getSharedPtr(), gctx, fourPosition,
          momentum.normalized(), a.charge() / momentum.norm(), a.covariance(),
          a.particleHypothesis());

      if (!res.ok()) {
        throw std::runtime_error("Invalid bound track parameters");
      }
      return res.value();
    } else {
      Acts::Vector4 fourPosition = a.fourPosition() + b.fourPosition();
      return TrackParameters(fourPosition, momentum.normalized(),
                             a.charge() / momentum.norm(), a.covariance(),
                             a.particleHypothesis());
    }
  }

  /// Grids to accumulate IP and reference
  /// layer track parameters
  LookupGrid m_grid;

  /// Mutex for protecting grid access
  std::mutex m_gridMutex;

  /// Map to keep the accumulation count
  /// in the occupied grid bins
  std::map<std::array<std::size_t, LookupGrid::DIM>, std::size_t> m_countGrid;
};

}  // namespace Acts
