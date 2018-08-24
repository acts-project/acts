// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// INavigationEngine.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"

namespace Acts {

class TrackingGeometry;

using ExCellCharged = ExtrapolationCell<TrackParameters>;
using ExCellNeutral = ExtrapolationCell<NeutralParameters>;

/// @class INavigationEngine
///
/// Extrapolation engine interface for Charged and Neutral parameters,
/// it serves as the Master extrapolation interface but also as the specialized
/// extrapolation engine ones.
///
class INavigationEngine
{
public:
  /// Virtual destructor
  virtual ~INavigationEngine() {}
  /// Resolve the boundary situation - for charged particles
  ///
  /// @param ecCell is the charged extrapolation cell
  /// @param dir is the additional direction prescription
  ///
  /// @return is a extrapolation code indication
  virtual ExtrapolationCode
  resolveBoundary(ExCellCharged&      ecCell,
                  NavigationDirection dir = forward) const = 0;

  /// Resolve the boundary situation - for neutral particles
  ///
  /// @param enCell is the neutral extrapolation cell
  /// @param dir is the additional direction prescription
  ///
  /// @return is a extrapolation code indication
  virtual ExtrapolationCode
  resolveBoundary(ExCellNeutral&      enCell,
                  NavigationDirection dir = forward) const = 0;

  /// Resolve the boundary situation - for charged particles
  ///
  /// @param ecCell is the charged extrapolation cell
  /// @param dir is the additional direction prescription
  /// @param noLoop is a loop protection
  /// @todo check with sharka where this is used
  ///
  /// @return is a extrapolation code indication
  virtual ExtrapolationCode
  resolvePosition(ExCellCharged&      ecCell,
                  NavigationDirection dir    = forward,
                  bool                noLoop = false) const = 0;

  /// Resolve the boundary situation - for neutral particles
  ///
  /// @param enCell is the neutral extrapolation cell
  /// @param dir is the additional direction prescription
  /// @param noLoop is a loop protection
  /// @todo check with sharka where this is used
  ///
  /// @return is a extrapolation code indication
  virtual ExtrapolationCode
  resolvePosition(ExCellNeutral&      enCell,
                  NavigationDirection dir    = forward,
                  bool                noLoop = false) const = 0;

protected:
  std::string m_sopPrefix;   ///< prefix for screen output
  std::string m_sopPostfix;  ///< prefix for screen output
};

}  // namespace