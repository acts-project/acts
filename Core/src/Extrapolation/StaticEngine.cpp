// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticEngine.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ACTS/Extrapolation/StaticEngine.hpp"

// constructor
Acts::StaticEngine::StaticEngine(const Acts::StaticEngine::Config& seConfig,
                                 std::unique_ptr<const Logger>     logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(seConfig);
}

// destructor
Acts::StaticEngine::~StaticEngine()
{
}

// configuration
void
Acts::StaticEngine::setConfiguration(const Acts::StaticEngine::Config& seConfig)
{
  // steering of the screen outoput (SOP)
  IExtrapolationEngine::m_sopPrefix  = seConfig.prefix;
  IExtrapolationEngine::m_sopPostfix = seConfig.postfix;
  // copy the configuration
  m_cfg = seConfig;
}

void
Acts::StaticEngine::setLogger(std::unique_ptr<const Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

// charged extrapolation /
Acts::ExtrapolationCode
Acts::StaticEngine::extrapolate(ExCellCharged&       ecCharged,
                                const Surface*       sf,
                                const BoundaryCheck& bcheck) const
{
  return extrapolateT<TrackParameters>(
      ecCharged, sf, ecCharged.propDirection, bcheck);
}

// neutral extrapolation
Acts::ExtrapolationCode
Acts::StaticEngine::extrapolate(ExCellNeutral&       ecNeutral,
                                const Surface*       sf,
                                const BoundaryCheck& bcheck) const
{
  return extrapolateT<NeutralParameters>(
      ecNeutral, sf, ecNeutral.propDirection, bcheck);
}
