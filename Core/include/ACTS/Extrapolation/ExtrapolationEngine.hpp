// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ExtrapolationEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATION_EXTRAPOLATIONENGINE_H
#define ACTS_EXTRAPOLATION_EXTRAPOLATIONENGINE_H 1

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/INavigationEngine.hpp"
#include "ACTS/Extrapolation/IPropagationEngine.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

class TrackingGeometry;
class Surface;
class BoundaryCheck;

/// @class ExtrapolationEngine
///
/// Master extrapolation engine for extrapolation through the TrackingGeometry.
///
/// It delegates the extrapolation to optimised engines, handing over the ExtrapolationCell
/// as internal cache.
///
/// There are identical interfaces for charged and neutral track parameters.
///
/// Providing a destination surface is optional, if no destination surface is given the extrapolation
/// process can be stopped by other directives, e.g. stopping at a certain path limit, material limit
/// or with a change of detector signature.
///
class ExtrapolationEngine : virtual public IExtrapolationEngine
{
public:
  /// @struct Config
  ///
  /// Configuration struct to be used for this ExtrapolationEngine.
  /// This has to be prepared by the framework/main and can be either given to
  ///the ExtrapolationEngine at construction or
  /// with the setConfig method.
  ///
  struct Config
  {
    /// the default logger
    std::shared_ptr<Logger> logger = getDefaultLogger("ExtrapolationEngine", Logging::INFO);
    /// the tracking geometry
    std::shared_ptr<const TrackingGeometry>   trackingGeometry = nullptr; 
    /// the list of extrapolation engines 
    std::vector<std::shared_ptr<IExtrapolationEngine>> extrapolationEngines;  
    /// the helper propagator for navigation initialization  
    std::shared_ptr<IPropagationEngine>       propagationEngine; 
    /// the navigation engine 
    std::shared_ptr<INavigationEngine>        navigationEngine;  
    std::string prefix;    ///< output prefix
    std::string postfix;   ///< output postfix
    std::string name;      ///< name of the tool
              
    Config()
      : logger(getDefaultLogger("ExtrapolationEngine", Logging::INFO))
      , trackingGeometry(nullptr)
      , extrapolationEngines()
      , propagationEngine(nullptr)
      , navigationEngine(nullptr)
      , prefix("[ME] - ")
      , postfix(" - ")
      , name("Anonymous")
    {
    }
  };

  /// Constructor 
  /// @param eeConfig is the configuration struct for this engine 
  ExtrapolationEngine(const Config& eeConfig);

  /// Destructor 
  ~ExtrapolationEngine();

  using IExtrapolationEngine::extrapolate;

  /// charged extrapolation - public interface 
  /// @param ecCharged is the charged extrapolation cell that holds the cache
  /// @param sf is the (optional) destinaton surface
  /// @param bchk is the boudnary check directive @TODO shift to cell after splitting
  ExtrapolationCode
  extrapolate(ExCellCharged&       ecCharged,
              const Surface*       sf     = nullptr,
              const BoundaryCheck& bchk = true) const final;

  /// neutral extrapolation - public interface 
  /// @param ecNeutral is the neutral extrapolation cell that holds the cache
  /// @param sf is the (optional) destinaton surface
  /// @param bchk is the boudnary check directive @TODO shift to cell after splitting
  ExtrapolationCode
  extrapolate(ExCellNeutral&       ecNeutral,
              const Surface*       sf     = nullptr,
              const BoundaryCheck& bchk = true) const final;

  /// define for which GeometrySignature this extrapolator is valid 
  ///  - this is GLOBAL
  GeometryType
  geometryType() const final;

  /// Set configuration method 
  void
  setConfiguration(const Config& eeConfig);

  /// Get configuration method 
  Config
  getConfiguration() const;

protected:
  /// ExtrapolationEngine config object 
  Config m_config;

private:
  /// Private access to the logging instance
  const Logger&
  logger() const
  {
    return *m_config.logger;
  }
  
  /// main extrapolation method, templated to chared/neutral 
  /// @paramt eCell ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param dir is the additional direction prescription
  /// @param bchk is the boudnary check directive @TODO shift to cell after splitting
  template <class T>
  ExtrapolationCode
  extrapolateT(ExtrapolationCell<T>& eCell,
               const Surface*        sf     = nullptr,
               PropDirection         dir    = alongMomentum,
               const BoundaryCheck&  bcheck = true) const;

  /// main extrapolation method, templated to chared/neutral 
  /// @paramt eCell ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param dir is the additional direction prescription
  template <class T>
  ExtrapolationCode
  initNavigation(ExtrapolationCell<T>& eCell,
                 const Surface*        sf  = nullptr,
                 PropDirection         dir = alongMomentum) const;
};

/// Return the geometry type, it's the master 
inline GeometryType
ExtrapolationEngine::geometryType() const
{
  return Acts::Master;
}

/// Return the configuration object 
inline ExtrapolationEngine::Config
ExtrapolationEngine::getConfiguration() const
{
  return m_config;
}

}  // end of namespace

#include "ACTS/Extrapolation/detail/ExtrapolationEngine.ipp"

#endif  // ACTS_EXTRAPOLATION_EXTRAPOLATIONENGINE_H
