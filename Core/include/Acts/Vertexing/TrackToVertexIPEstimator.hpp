// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

struct ImpactParametersAndSigma {
  double IPd0 = 0.;
  double IPz0 = 0.;
  double IPz0SinTheta = 0.;
  double sigmad0 = 0.;
  double sigmaz0 = 0.;
  double sigmaz0SinTheta = 0.;
  double PVsigmad0 = 0.;
  double PVsigmaz0 = 0.;
  double PVsigmaz0SinTheta = 0.;
};

template <typename input_track_t, typename propagator_t,
          typename action_list_t = ActionList<>,
          typename aborter_list_t = AbortList<>>

/// @class TrackToVertexIPEstimator estimates the impact parameters and their
/// errors of a given track w.r.t. a vertex
class TrackToVertexIPEstimator {
 public:
  /// @brief Default constructor
  ///
  /// @param logger Logging instance
  TrackToVertexIPEstimator(std::unique_ptr<const Logger> logger =
                               getDefaultLogger("TrackToVertexIPEstimator",
                                                Logging::INFO))
      : m_logger(std::move(logger)) {}

  /// @brief Move constructor
  TrackToVertexIPEstimator(TrackToVertexIPEstimator&& other)
      : m_logger(std::move(other.m_logger)) {}

  /// @brief Estimates the impact parameters and their errors of a given
  /// track w.r.t. a vertex by propagating the trajectory state
  /// towards the vertex position.
  ///
  /// @param gctx The Geometry Context
  /// @param mctx The MagneticField Context
  /// @param track Track to estimate IP from
  /// @param vtx Vertex the track belongs to
  /// @param propagator Propagator
  Result<std::unique_ptr<ImpactParametersAndSigma>> estimate(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const BoundParameters& track, const Vertex<input_track_t>& vtx,
      const propagator_t& propagator) const;

 private:
  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "TrackToVertexIPEstimator.ipp"