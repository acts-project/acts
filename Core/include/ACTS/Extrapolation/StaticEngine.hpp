// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATION_STATICENGINE_H
#define ACTS_EXTRAPOLATION_STATICENGINE_H 1

#ifndef ACTS_EXTRAPOLATION_OUTPUTHELPER
#define ACTS_EXTRAPOLATION_OUTPUTHELPER 1
#define OH_CHECKFOUND(object) (object ? "found" : "not found")
#endif

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/IMaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/INavigationEngine.hpp"
#include "ACTS/Extrapolation/IPropagationEngine.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

/// @class StaticEngine
///
/// Extrapolation engine for static layer & volume setup.
///
/// This engine relies on the fact that every position in a static layer setup
/// can
/// be uniquely associated to a layer (NavigationLayer or real physical layer),
/// and thus to a voume at navigation level. The extrapolation process within a
/// fully static setup is then realised as a step from layer to layer within a
/// volume,
/// and from volume to volume at a higher level.
///
/// The entire code is written as a template in either charged or neutral
/// parameters.
///

class StaticEngine : virtual public IExtrapolationEngine
{
public:
  /// @enum ResolveLayerType
  /// - use for code readability
  ///
  enum ResolveLayerType {
    StartLayer               = 0,
    NavigationLayer          = 1,
    PassThroughLayer         = 2,
    SubStructureLayer        = 3,
    DestinationLayer         = 4,
    StartAndDestinationLayer = 6,
    UndefinedLayer           = 5
  };

  /// @struct Config
  /// Configuration struct of the StaticEngine,
  /// holds the helper engines that can be configured
  struct Config
  {
    /// the used propagation engine
    std::shared_ptr<const IPropagationEngine> propagationEngine = nullptr;
    /// the navigation engine to resolve the boundary
    std::shared_ptr<const INavigationEngine> navigationEngine = nullptr;
    /// the material effects updated
    std::shared_ptr<const IMaterialEffectsEngine> materialEffectsEngine
        = nullptr;
    /// output prefix
    std::string prefix = "[SE] - ";
    /// output postfix
    std::string postfix = " - ";
  };

  /// Constructor
  ///
  /// @param seConfig is the configuration struct
  /// @param logger logging instance
  StaticEngine(const Config&           seConfig,
               std::unique_ptr<Logger> logger
               = getDefaultLogger("StaticEngine", Logging::INFO));

  /// Destructor
  ~StaticEngine();

  using IExtrapolationEngine::extrapolate;

  /// Main extrapolation method, templated to chared
  ///
  /// @param ecCharged ist he extrapolation cell for charged parameters
  /// @param sf is the (optional) destinaton surface
  /// @param bcheck is the boudnary check directive @todo shift to cell after
  /// splitting
  ///
  /// @return is a extrapolation code indication
  ExtrapolationCode
  extrapolate(ExCellCharged&       ecCharged,
              const Surface*       sf     = 0,
              const BoundaryCheck& bcheck = true) const final;

  /// Main extrapolation method, templated to neutral
  ///
  /// @param ecNeutral ist he extrapolation cell for neutral parameters
  /// @param sf is the (optional) destinaton surface
  /// @param bcheck is the boudnary check directive @todo shift to cell after
  /// splitting
  ///
  /// @return is a extrapolation code indication
  ExtrapolationCode
  extrapolate(ExCellNeutral&       ecNeutral,
              const Surface*       sf     = 0,
              const BoundaryCheck& bcheck = true) const final;

  /// Define for which geometry type this extrapolator is valid
  ///
  /// @return this retursn static for this engine
  GeometryType
  geometryType() const final;

  /// Set configuration method
  ///
  /// @param seConfig the configuration object to be set
  void
  setConfiguration(const Config& seConfig);

  /// Get configuration method
  Config
  getConfiguration() const;

  /// Set logging instance
  ///
  /// @param logger is the logging instance to be set
  void
  setLogger(std::unique_ptr<Logger> logger);

protected:
  /// Configuration struct
  Config m_cfg;

private:
  /// Private access to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  std::unique_ptr<Logger> m_logger;

  /// Main loop extrapolation method
  ///
  /// @param eCell ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param dir is the additional direction prescription
  /// @param bcheck is the boudnary check directive @todo shift to cell after
  /// splitting
  ///
  /// @return is a extrapolation code indication
  template <class T>
  ExtrapolationCode
  extrapolateT(ExtrapolationCell<T>& eCell,
               const Surface*        sf     = 0,
               PropDirection         dir    = alongMomentum,
               const BoundaryCheck&  bcheck = true) const;

  /// Init Navigation for static setup
  ///
  /// @param eCell ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param dir is the additional direction prescription
  /// @param bcheck is the boudnary check directive @todo shift to cell after
  /// splitting
  ///
  /// @return is a extrapolation code indication
  template <class T>
  ExtrapolationCode
  initNavigationT(ExtrapolationCell<T>& eCell,
                  const Surface*        sf     = 0,
                  PropDirection         dir    = alongMomentum,
                  const BoundaryCheck&  bcheck = true) const;

  /// Main static layer handling
  ///
  /// @param eCell ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param dir is the additional direction prescription
  /// @param bcheck is the boudnary check directive @todo shift to cell after
  /// splitting
  ///
  /// @return is a extrapolation code indication
  template <class T>
  ExtrapolationCode
  handleLayerT(ExtrapolationCell<T>& eCell,
               const Surface*        sf     = 0,
               PropDirection         dir    = alongMomentum,
               const BoundaryCheck&  bcheck = true) const;

  /// Main sub structure layer handling
  ///
  /// @param eCell ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param dir is the additional direction prescription
  /// @param bcheck is the boudnary check directive @todo shift to cell after
  /// splitting
  /// @param hasSubStructure is an indicator whether the layer has sub structure
  /// which needs to be resolved
  /// @param isStartLayer is and indicator whether the layer
  /// is the start layer
  /// @param isDestinationLayer is and indicator whether the layer
  /// is the destination layer
  ///
  /// @return is a extrapolation code indication
  template <class T>
  ExtrapolationCode
  resolveLayerT(ExtrapolationCell<T>& eCell,
                const Acts::Surface*  sf,
                PropDirection         dir                = alongMomentum,
                const BoundaryCheck&  bcheck             = true,
                bool                  hasSubStructure    = false,
                bool                  isStartLayer       = false,
                bool                  isDestinationLayer = false) const;

  /// Handle the failure - as configured
  ///
  /// @param eCode is the extrapolation code at entry
  /// @param eCell ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param dir is the additional direction prescription
  /// @param bcheck is the boudnary check directive @todo shift to cell after
  /// splitting
  ///
  /// @return is a extrapolation code indication at exit
  template <class T>
  ExtrapolationCode
  handleReturnT(ExtrapolationCode     eCode,
                ExtrapolationCell<T>& eCell,
                const Surface*        sf     = 0,
                PropDirection         dir    = alongMomentum,
                const BoundaryCheck&  bcheck = true) const;
};

/// Return the geomtry type
inline GeometryType
StaticEngine::geometryType() const
{
  return Acts::Static;
}

/// Return the configuration object
inline StaticEngine::Config
StaticEngine::getConfiguration() const
{
  return m_cfg;
}

}  // end of namespace

#include "ACTS/Extrapolation/detail/StaticEngine.ipp"

#endif  // ACTS_EXTRAPOLATION_STATICENGINE_H
