// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialMapper.hpp, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class TrackingGeometry;

/// @brief SurfaceMaterialMapper
///
/// This is the main feature tool to map material information
/// from a 3D geometry onto the TrackingGeometry with it's surface
/// material description.
///
/// The process runs as such:
///
///  1) TrackingGeometry is parsed and for each Surface with
///     with SurfaceMaterialProxy a local store is initialized
///     the identification is done hereby through the Surface::GeometryID
///
///  2) A Cache is generated that is used to keep the filling thread local,\
  ///     the filling is protected with std::mutex
///
///  3) A number of N material tracks is read in, each track has :
///       origin, direction, material steps < position, step length, x0, l0, a,
///       z, rho >
///
///       for each track:
///          surfaces along the origin/direction path are collected
///          the closest material steps are assigned
///
///  4) Each 'hit' bin per event is counted and averaged at the end of the run
///
class SurfaceMaterialMapper
{
public:
  using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;

  /// @struct Config
  ///
  /// Nested Configuration struct for the material mapper
  struct Config
  {
    std::array<double, 2> etaRange = {{-6., 6.}};
  };

  /// Delete the Default constructor
  SurfaceMaterialMapper() = delete;

  /// Constructor with config object
  SurfaceMaterialMapper(const Config&                 cfg,
                        StraightLinePropagator        propagator,
                        std::unique_ptr<const Logger> log
                        = getDefaultLogger("SurfaceMaterialMapper",
                                           Logging::INFO));

  /// Default destructor
  ~SurfaceMaterialMapper() = default;

  /// Process a single 'track'
  void
  processTrack() const;

private:
  /// Standard logger method
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// The configuration object
  Config m_cfg;

  /// The logging instance
  std::unique_ptr<const Logger> m_logger;

  /// The straight line propagator
  StraightLinePropagator m_propagator;
};
}
