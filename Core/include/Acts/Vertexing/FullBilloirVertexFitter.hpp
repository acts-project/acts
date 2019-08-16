// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizerConcept.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitterOptions.hpp"
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
/// @tparam bfield_t Magnetic field type
/// @tparam input_track_t Track object type
/// @tparam linearizer_t Track linearizer type
template <typename bfield_t, typename input_track_t,
          typename linearizer_t = HelicalTrackLinearizer<bfield_t>>
class FullBilloirVertexFitter {
  static_assert(LinearizerConcept<linearizer_t>,
                "Linearizer does not fulfill linearizer concept.");

 public:
  using InputTrack_t = input_track_t;
  using Propagator_t = typename linearizer_t::Propagator_t;
  using Linearizer_t = linearizer_t;
  using Bfield_t = bfield_t;

  struct Config {
    /// Magnetic field
    bfield_t bField;

    /// Maximum number of interations in fitter
    int maxIterations = 5;

    /// Set up factory for linearizing tracks
    typename linearizer_t::Config linConfig;
    linearizer_t linearizer;

    /// Propagator
    Propagator_t propagator;

    /// Constructor with propagator input
    Config(const bfield_t& bIn, const Propagator_t& propagatorIn)
        : bField(bIn),
          linConfig(bIn),
          linearizer(linConfig),
          propagator(propagatorIn) {}

    /// Constructor with default propagator
    template <typename T = Propagator_t,
              std::enable_if_t<
                  std::is_same<T, Propagator<EigenStepper<bfield_t>>>::value,
                  int> = 0>
    Config(const bfield_t& bIn)
        : bField(bIn),
          linConfig(bIn),
          linearizer(linConfig),
          propagator(
              Propagator<EigenStepper<bfield_t>>(EigenStepper<bfield_t>(bIn))) {
    }
  };

  /// @brief Constructor used if input_track_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  template <typename T = input_track_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  FullBilloirVertexFitter(const Config& cfg)
      : m_cfg(cfg), extractParameters([](T params) { return params; }) {}

  /// @brief Constructor for user-defined input_track_t type =! BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from input_track_t object
  FullBilloirVertexFitter(const Config& cfg,
                          std::function<BoundParameters(input_track_t)> func)
      : m_cfg(cfg), extractParameters(func) {}

  /// @brief Fit method, fitting vertex for provided tracks with constraint
  ///
  /// @param paramVector Vector of track objects to fit vertex to
  /// @param vFitterOptions Vertex fitter options
  ///
  /// @return Fitted vertex
  Result<Vertex<input_track_t>> fit(
      const std::vector<input_track_t>& paramVector,
      const VertexFitterOptions<input_track_t>& vFitterOptions) const;

 private:
  /// Configuration object
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// input_track_t objects are BoundParameters by default, function to be
  /// overwritten to return BoundParameters for other input_track_t objects.
  ///
  /// @param params input_track_t object to extract track parameters from
  std::function<BoundParameters(input_track_t)> extractParameters;

  /// @brief Function to correct 2-pi periodicity for phi and theta
  ///
  /// @param phiIn Phi
  /// @param thetaIn Theta
  ///
  /// @return Pair of (corrected phi, corrected theta)
  std::pair<double, double> correctPhiThetaPeriodicity(double phiIn,
                                                       double thetaIn) const;
};

}  // namespace Acts

#include "FullBilloirVertexFitter.ipp"
