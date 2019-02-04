// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 Acts project team
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

template <typename BField, typename InputTrack>
class FullBilloirVertexFitter : public IVertexFitter<InputTrack>
{
public:
  struct Config
  {
    /// Magnetic field
    BField bField;
    /// Starting point of vertex fit
    Vector3D startingPoint;
    /// Maximum number of interations in fitter
    int maxIterations;

    /// Constructor with default number of iterations and starting point
    Config(BField bIn) : bField(std::move(bIn)), maxIterations(5) {}
  };

  /// Constructor with explicit config
  FullBilloirVertexFitter(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Fit method, fitting vertex for provided tracks with constraint
  ///
  /// @param paramVector Vector of track objects to fit vertex to
  /// @param startingPoint Constraint of the fit, position of constraint is
  /// starting point
  ///
  /// @return Fitted vertex
  Vertex<InputTrack>
  fit(const std::vector<InputTrack>& paramVector,
      Vertex<InputTrack>             constraint = Vertex<InputTrack>(Vector3D(0.,0.,0.))) const;

private:
  /// Configuration object
  Config m_cfg;
};

}  // namespace Acts

#include "FullBilloirVertexFitter.ipp"