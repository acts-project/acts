// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// INavigationEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONINTERFACES_INAVIGATIONENGINE_H
#define ACTS_EXTRAPOLATIONINTERFACES_INAVIGATIONENGINE_H

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"

namespace Acts {

class TrackingGeometry;

typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

/** @class INavigationEngine

    Extrapolation engine interface for Charged and Neutral parameters,
    it serves as the Master extrapolation interface but also as the specialized
    extrapolation engine ones.

*/

class INavigationEngine
{
public:
  /** Virtual destructor */
  virtual ~INavigationEngine() {}
  /** resolve the boundary situation - for charged particles */
  virtual ExtrapolationCode
  resolveBoundary(ExCellCharged& ecCell,
                  PropDirection  dir = alongMomentum) const = 0;

  /** resolve the boundary situation - for neutral particles */
  virtual ExtrapolationCode
  resolveBoundary(ExCellNeutral& enCell,
                  PropDirection  dir = alongMomentum) const = 0;

  /** resolve the position - for charged particles */
  virtual ExtrapolationCode
  resolvePosition(ExCellCharged& ecCell,
                  PropDirection  dir    = alongMomentum,
                  bool           noLoop = false) const = 0;

  /** resolve the position - for neutral particles */
  virtual ExtrapolationCode
  resolvePosition(ExCellNeutral& enCell,
                  PropDirection  dir    = alongMomentum,
                  bool           noLoop = false) const = 0;

protected:
  //!< SCREEN output formatting  (SOP) - unify amongst extrapolation engines
  std::string m_sopPrefix;   //!< prefix for screen output
  std::string m_sopPostfix;  //!< prefix for screen output
};

}  // end of namespace

#endif  // ACTS_EXTRAPOLATIONINTERFACES_INAVIGATIONENGINE_H
