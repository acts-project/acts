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

#ifndef ACTS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H
#define ACTS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H 1

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/IMaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/INavigationEngine.hpp"
#include "ACTS/Extrapolation/IPropagationEngine.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Volumes/BoundarySurface.hpp"

namespace Acts {

class TrackingGeometry;

/** @class StaticNavigationEntine

    The static navigation engine for finding the next volume,
    propagate to the boundary, can be shared with other engines that have a
   static frame.

  */
class StaticNavigationEngine : virtual public INavigationEngine
{
public:
  /** @struct Config
      configuration struct for the StaticNavigationEngine
  */
  struct Config
  {
    std::shared_ptr<Logger> logger;

    std::shared_ptr<IPropagationEngine>
        propagationEngine;  //!< the used propagation engine
    std::shared_ptr<IMaterialEffectsEngine>
        materialEffectsEngine;  //!< the material effects updated
    std::shared_ptr<const TrackingGeometry>
        trackingGeometry;  //!< the tracking geometry used by the navigator

    std::string prefix;   //!< output prefix
    std::string postfix;  //!< output postfix
    std::string name;     //!< name of this engine

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

  /** Constructor */
  StaticNavigationEngine(const Config& snConfig);

  /** Destructor */
  ~StaticNavigationEngine();

  /** avoid method shaddowing */
  using INavigationEngine::resolveBoundary;
  using INavigationEngine::resolvePosition;

  /** resolve the boundary situation - for charged particles */
  ExtrapolationCode
  resolveBoundary(ExCellCharged& eCell,
                  PropDirection  dir = alongMomentum) const final;

  /** resolve the boundary situation - for neutral particles */
  ExtrapolationCode
  resolveBoundary(ExCellNeutral& eCelll,
                  PropDirection  dir = alongMomentum) const final;

  /** resolve the boundary situation - for charged particles */
  ExtrapolationCode
  resolvePosition(ExCellCharged& eCell,
                  PropDirection  dir    = alongMomentum,
                  bool           noLoop = false) const final;

  /** resolve the boundary situation - for neutral particles */
  ExtrapolationCode
  resolvePosition(ExCellNeutral& eCelll,
                  PropDirection  dir    = alongMomentum,
                  bool           noLoop = false) const final;

  /** Set configuration method */
  void
  setConfiguration(const Config& meConfig);

  /** Get configuration method */
  Config
  getConfiguration() const;

protected:
  /** the configuration member of the static navigation engine */
  Config m_config;

private:
  const Logger&
  logger() const
  {
    return *m_config.logger;
  }
  /** resolve the boundary situation */
  template <class T>
  ExtrapolationCode
  resolveBoundaryT(ExtrapolationCell<T>& eCell,
                   PropDirection         dir = alongMomentum) const;

  /** resolve position */
  template <class T>
  ExtrapolationCode
  resolvePositionT(ExtrapolationCell<T>& eCell,
                   PropDirection         dir    = alongMomentum,
                   bool                  noLoop = false) const;

  /** deal with the boundary Surface - called by resolveBoundary */
  template <class T>
  ExtrapolationCode
  handleBoundaryT(ExtrapolationCell<T>&                  eCell,
                  const BoundarySurface<TrackingVolume>& bSurfaceTV,
                  PropDirection                          dir = alongMomentum,
                  bool                                   stepout = false) const;
};

/** Return the configuration object */
inline StaticNavigationEngine::Config
StaticNavigationEngine::getConfiguration() const
{
  return m_config;
}

}  // end of namespace

//!< define the templated function
#include "ACTS/Extrapolation/detail/StaticNavigationEngine.ipp"

#endif  // ACTS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H
