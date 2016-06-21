// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATION_STATICNAVIGATIONENGINE_H
#define ACTS_EXTRAPOLATION_STATICNAVIGATIONENGINE_H 1

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/IMaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/INavigationEngine.hpp"
#include "ACTS/Extrapolation/IPropagationEngine.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Volumes/BoundarySurfaceT.hpp"

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
  /// NEsted configuration struct for the StaticNavigationEngine
  struct Config
  {
    std::shared_ptr<Logger> logger;

    /// the used propagation engine
    std::shared_ptr<IPropagationEngine>       propagationEngine;  
    /// the material effects updator
    std::shared_ptr<IMaterialEffectsEngine>   materialEffectsEngine;  
    /// the tracking geometry cache
    std::shared_ptr<const TrackingGeometry>   trackingGeometry;  

    std::string prefix;   ///< output prefix
    std::string postfix;  ///< output postfix
    std::string name;     ///< name of this engine

    Config()
      : logger(getDefaultLogger("StaticNavigationEngine", Logging::INFO))
      , propagationEngine(nullptr)
      , materialEffectsEngine(nullptr)
      , trackingGeometry(nullptr)
      , prefix("[SN] - ")
      , postfix(" - ")
      , name("Anonymous")
    {
    }
  };

  /// Constructor 
  /// @param snConfig is the configuration struct to steer behaviour
  StaticNavigationEngine(const Config& snConfig);

  /// Destructor 
  ~StaticNavigationEngine();

  /// avoid method shaddowing 
  using INavigationEngine::resolveBoundary;
  using INavigationEngine::resolvePosition;

  /// resolve the boundary situation - for charged particles 
  /// @param ecCell is the charged extrapolation cell
  /// @param dir is the additional direction prescription
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolveBoundary(ExCellCharged& eCell,
                  PropDirection  dir = alongMomentum) const final;

  /// resolve the boundary situation - for neutral particles 
  /// @param ecCell is the neutral extrapolation cell
  /// @param dir is the additional direction prescription
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolveBoundary(ExCellNeutral& eCelll,
                  PropDirection  dir = alongMomentum) const final;

  /// resolve the boundary situation - for charged particles 
  /// @param ecCell is the charged extrapolation cell
  /// @param dir is the additional direction prescription
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolvePosition(ExCellCharged& eCell,
                  PropDirection  dir    = alongMomentum,
                  bool           noLoop = false) const final;

  /// resolve the boundary situation - for neutral particles 
  /// @param ecCell is the neutral extrapolation cell
  /// @param dir is the additional direction prescription
  /// @return is a extrapolation code indication
  ExtrapolationCode
  resolvePosition(ExCellNeutral& eCelll,
                  PropDirection  dir    = alongMomentum,
                  bool           noLoop = false) const final;

  /// Set configuration method 
  void
  setConfiguration(const Config& meConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

protected:
  /// the configuration member of the static navigation engine 
  Config m_config;

private:
  /// Private access to the logging instance
  const Logger&
  logger() const
  {
    return *m_config.logger;
  }

  /// resolve the boundary situation 
  template <class T>
  ExtrapolationCode
  resolveBoundaryT(ExtrapolationCell<T>& eCell,
                   PropDirection         dir = alongMomentum) const;

  /// resolve position 
  template <class T>
  ExtrapolationCode
  resolvePositionT(ExtrapolationCell<T>& eCell,
                   PropDirection         dir    = alongMomentum,
                   bool                  noLoop = false) const;

  /// deal with the boundary Surface - called by resolveBoundary 
  template <class T>
  ExtrapolationCode
  handleBoundaryT(ExtrapolationCell<T>&                  eCell,
                  const BoundarySurfaceT<TrackingVolume>& bSurfaceTV,
                  PropDirection                          dir = alongMomentum,
                  bool                                   stepout = false) const;
};

inline StaticNavigationEngine::Config
StaticNavigationEngine::getConfiguration() const
{
  return m_config;
}

}  // end of namespace

#include "ACTS/Extrapolation/detail/StaticNavigationEngine.ipp"

#endif  // ACTS_EXTRAPOLATION_STATICNAVIGATIONENGINE_H
