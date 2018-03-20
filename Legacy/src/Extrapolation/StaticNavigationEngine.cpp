// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.cpp, Acts project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "Acts/Extrapolation/StaticNavigationEngine.hpp"

// constructor
Acts::StaticNavigationEngine::StaticNavigationEngine(
    const Acts::StaticNavigationEngine::Config& snConfig,
    std::unique_ptr<const Logger>               logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(snConfig);
}

// destructor
Acts::StaticNavigationEngine::~StaticNavigationEngine()
{
}

// configuration
void
Acts::StaticNavigationEngine::setConfiguration(
    const Acts::StaticNavigationEngine::Config& snConfig)
{
  // steering of the screen outoput (SOP)
  INavigationEngine::m_sopPrefix  = snConfig.prefix;
  INavigationEngine::m_sopPostfix = snConfig.postfix;
  // copy the configuration
  m_cfg = snConfig;
}

void
Acts::StaticNavigationEngine::setLogger(std::unique_ptr<const Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

// charged situation
Acts::ExtrapolationCode
Acts::StaticNavigationEngine::resolveBoundary(Acts::ExCellCharged& ecCharged,
                                              PropDirection        dir) const
{
  return resolveBoundaryT<Acts::TrackParameters>(ecCharged, dir);
}

// neutral situation
Acts::ExtrapolationCode
Acts::StaticNavigationEngine::resolveBoundary(Acts::ExCellNeutral& ecNeutral,
                                              PropDirection        dir) const
{
  return resolveBoundaryT<Acts::NeutralParameters>(ecNeutral, dir);
}

// charged situation
Acts::ExtrapolationCode
Acts::StaticNavigationEngine::resolvePosition(Acts::ExCellCharged& ecCharged,
                                              PropDirection        dir,
                                              bool                 noLoop) const
{
  return resolvePositionT<Acts::TrackParameters>(ecCharged, dir, noLoop);
}

// neutral situation
Acts::ExtrapolationCode
Acts::StaticNavigationEngine::resolvePosition(Acts::ExCellNeutral& ecNeutral,
                                              PropDirection        dir,
                                              bool                 noLoop) const
{
  return resolvePositionT<Acts::NeutralParameters>(ecNeutral, dir, noLoop);
}
