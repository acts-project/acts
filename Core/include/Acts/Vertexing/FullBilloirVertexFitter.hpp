// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizerConcept.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
namespace Acts {

/// @class FullBilloirVertexFitter
///
/// @brief Vertex fitter class implementing the Billoir vertex fitter
///
/// This class implements the Billoir vertex fitter:
///
/// Fast vertex fitting with a local parametrization of tracks
/// Author(s) Billoir, P ; Qian, S
/// In: Nucl. Instrum. Methods Phys. Res., A 311 (1992) 139-150
/// DOI 10.1016/0168-9002(92)90859-3
///
/// @tparam input_track_t Track object type
/// @tparam linearizer_t Track linearizer type
template <typename input_track_t, typename linearizer_t>
class FullBilloirVertexFitter {
  static_assert(LinearizerConcept<linearizer_t>,
                "Linearizer does not fulfill linearizer concept.");

 public:
  using InputTrack_t = input_track_t;
  using Propagator_t = typename linearizer_t::Propagator_t;
  using Linearizer_t = linearizer_t;

  struct State {
    /// @brief The state constructor
    ///
    /// @param fieldCache The magnetic field cache
    State(MagneticFieldProvider::Cache fieldCache)
        : linearizerState(std::move(fieldCache)) {}
    /// The linearizer state
    typename Linearizer_t::State linearizerState;
  };

  struct Config {
    /// Maximum number of interations in fitter
    int maxIterations = 5;
  };

  /// @brief Constructor used if input_track_t type == BoundTrackParameters
  ///
  /// @param cfg Configuration object
  template <
      typename T = input_track_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>
  FullBilloirVertexFitter(const Config& cfg)
      : m_cfg(cfg), extractParameters([](T params) { return params; }) {}

  /// @brief Constructor for user-defined input_track_t type =!
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from input_track_t
  /// object
  FullBilloirVertexFitter(
      const Config& cfg,
      std::function<BoundTrackParameters(input_track_t)> func)
      : m_cfg(cfg), extractParameters(func) {}

  /// @brief Fit method, fitting vertex for provided tracks with constraint
  ///
  /// @param paramVector Vector of track objects to fit vertex to
  /// @param linearizer The track linearizer
  /// @param vertexingOptions Vertexing options
  /// @param state The state object
  ///
  /// @return Fitted vertex
  Result<Vertex<input_track_t>> fit(
      const std::vector<const input_track_t*>& paramVector,
      const linearizer_t& linearizer,
      const VertexingOptions<input_track_t>& vertexingOptions,
      State& state) const;

 private:
  /// Configuration object
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// input_track_t objects are BoundTrackParameters by default, function to be
  /// overwritten to return BoundTrackParameters for other input_track_t
  /// objects.
  std::function<BoundTrackParameters(input_track_t)> extractParameters;
};

}  // namespace Acts

#include "FullBilloirVertexFitter.ipp"
