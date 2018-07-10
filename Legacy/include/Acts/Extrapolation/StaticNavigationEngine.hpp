// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"
#include "Acts/Extrapolation/IMaterialEffectsEngine.hpp"
#include "Acts/Extrapolation/INavigationEngine.hpp"
#include "Acts/Extrapolation/IPropagationEngine.hpp"
#include "Acts/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Volumes/BoundarySurfaceT.hpp"

namespace Acts {

class TrackingGeometry;

/// @class StaticNavigationEntine
///
/// The static navigation engine for finding the next volume,
/// propagate to the boundary, can be shared with other engines that have a
/// static frame.
///
class StaticNavigationEngine : virtual public INavigationEngine
{
public:
  /// @struct Config
  /// Nested configuration struct for the StaticNavigationEngine
  struct Config
  {
    /// the used propagation engine
    std::shared_ptr<const IPropagationEngine> propagationEngine = nullptr;
    /// the material effects updator
    std::shared_ptr<const IMaterialEffectsEngine> materialEffectsEngine
        = nullptr;
    /// the tracking geometry cache
    std::shared_ptr<const TrackingGeometry> trackingGeometry = nullptr;
    /// output prefix
    std::string prefix = "[SN] - ";
    /// output postfix
    std::string postfix = " - ";
  };

  /// Constructor
  ///
  /// @param snConfig is the configuration struct to steer behaviour
  /// @param logger logging instance
  StaticNavigationEngine(const Config&                 snConfig,
                         std::unique_ptr<const Logger> logger
                         = getDefaultLogger("StaticNavigationEngine",
                                            Logging::INFO));

  /// Destructor
  ~StaticNavigationEngine();

  /// avoid method shaddowing
  using INavigationEngine::resolveBoundary;
  using INavigationEngine::resolvePosition;

  /// Resolve the boundary situation - for charged particles
  ///
  /// @param ecCell is the charged extrapolation cell
  /// @param dir is the additional direction prescription
  ///
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolveBoundary(ExCellCharged&      ecCell,
                  NavigationDirection dir = forward) const final;

  /// Resolve the boundary situation - for neutral particles
  ///
  /// @param enCell is the neutral extrapolation cell
  /// @param dir is the additional direction prescription
  ///
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolveBoundary(ExCellNeutral&      enCell,
                  NavigationDirection dir = forward) const final;

  /// Resolve the boundary situation - for charged particles
  ///
  /// @param ecCell is the charged extrapolation cell
  /// @param dir is the additional direction prescription
  /// @param noLoop is a loop protection @todo check with ST
  ///
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolvePosition(ExCellCharged&      ecCell,
                  NavigationDirection dir    = forward,
                  bool                noLoop = false) const final;

  /// Resolve the boundary situation - for neutral particles
  ///
  /// @param enCell is the neutral extrapolation cell
  /// @param dir is the additional direction prescription
  /// @param noLoop is a loop protection @todo check with ST
  ///
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolvePosition(ExCellNeutral&      enCell,
                  NavigationDirection dir    = forward,
                  bool                noLoop = false) const final;

  /// Set configuration method
  ///
  /// @param snConfig the configuration object to be set
  void
  setConfiguration(const Config& snConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// Set logging instance
  ///
  /// @param logger the logging instance to be seet
  void
  setLogger(std::unique_ptr<const Logger> logger);

protected:
  /// the configuration member of the static navigation engine
  Config m_cfg;

private:
  /// Private access to the logging instance
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  std::unique_ptr<const Logger> m_logger;

  /// Resolve the boundary situation
  ///
  /// @param eCell the extrapolation
  /// @param dir the propagation direction
  ///
  /// @return is a extrapolation code indication
  template <class T>
  ExtrapolationCode
  resolveBoundaryT(ExtrapolationCell<T>& eCell,
                   NavigationDirection   dir = forward) const;

  /// Resolve position
  ///
  /// @param eCell the extrapolation
  /// @param dir the propagation direction
  /// @param noLoop
  /// @todo check with sharka
  ///
  /// @return is a extrapolation code indication
  template <class T>
  ExtrapolationCode
  resolvePositionT(ExtrapolationCell<T>& eCell,
                   NavigationDirection   dir    = forward,
                   bool                  noLoop = false) const;

  /// Deal with the boundary Surface - called by resolveBoundary
  ///
  /// @param eCell the extrapolation
  /// @param bSurfaceTV the boundary surface
  /// @param dir the propagation direction
  /// @param stepout  is a prescription to step out the volume
  ///
  /// @return is a extrapolation code indication
  template <class T>
  ExtrapolationCode
  handleBoundaryT(ExtrapolationCell<T>&                   eCell,
                  const BoundarySurfaceT<TrackingVolume>& bSurfaceTV,
                  NavigationDirection                     dir = forward,
                  bool stepout                                = false) const;
};

inline StaticNavigationEngine::Config
StaticNavigationEngine::getConfiguration() const
{
  return m_cfg;
}

}  // namespace

#include "Acts/Extrapolation/detail/StaticNavigationEngine.ipp"