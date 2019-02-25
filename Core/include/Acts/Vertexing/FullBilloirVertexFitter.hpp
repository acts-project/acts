// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Vertexing/IVertexFitter.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @class FullBilloirVertexFitter
///
/// @brief Vertex fitter class implementing the Billoir vertex fitter
///
/// This class implements the Billoir vertex fitter:
///
/// Fast vertex fitting with a local parametrization of tracks
/// Author(s)	Billoir, P ; Qian, S
/// In:	Nucl. Instrum. Methods Phys. Res., A 311 (1992) 139-150
/// DOI	10.1016/0168-9002(92)90859-3

template <typename BField, typename InputTrack, typename Propagator_t>
class FullBilloirVertexFitter : public IVertexFitter<InputTrack, Propagator_t>
{
public:
  struct Config
  {
    /// Magnetic field
    BField bField;

    /// Maximum number of interations in fitter
    int maxIterations = 5;

    /// Constructor with default number of iterations and starting point
    Config(BField bIn) : bField(std::move(bIn)) {}
  };

  /// @brief Constructor used if InputTrack type == BoundParameters
  ///
  /// @param cfg Configuration object
  template <typename T = InputTrack,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  FullBilloirVertexFitter(const Config& cfg)
    : m_cfg(cfg), extractParameters([&](T params) { return params; })
  {
  }

  /// @brief Constructor for user-defined InputTrack type =! BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from InputTrack object
  FullBilloirVertexFitter(const Config&                              cfg,
                          std::function<BoundParameters(InputTrack)> func)
    : m_cfg(cfg), extractParameters(func)
  {
  }

  /// @brief Fit method, fitting vertex for provided tracks with constraint
  ///
  /// @param paramVector Vector of track objects to fit vertex to
  /// @param propagator Propagator
  /// @param constraint Constraint of the fit, position of constraint is
  /// starting point
  ///
  /// @return Fitted vertex
  Vertex<InputTrack>
  fit(const std::vector<InputTrack>& paramVector,
      const Propagator_t&            propagator,
      Vertex<InputTrack>             constraint
      = Vertex<InputTrack>(Vector3D(0., 0., 0.))) const override;

private:
  /// Configuration object
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack objects are BoundParameters by default, function to be
  /// overwritten to return BoundParameters for other InputTrack objects.
  ///
  /// @param params InputTrack object to extract track parameters from
  std::function<BoundParameters(InputTrack)> extractParameters;
};

}  // namespace Acts

#include "FullBilloirVertexFitter.ipp"