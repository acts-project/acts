// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IExtrapolationEngine.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"
#include "Acts/Utilities/GeometrySignature.hpp"

namespace Acts {

using ExCellCharged = ExtrapolationCell<TrackParameters>;
using ExCellNeutral = ExtrapolationCell<NeutralParameters>;

class Surface;
class BoundaryCheck;

/// @class IExtrapolationEngine
///
/// Extrapolation engine interface for Charged and Neutral parameters,
/// it serves as the Master extrapolation interface but also as the specialized
/// extrapolation engine ones.
///
/// The ExtrapolationEngine is desinged as a thread safe const-correct service,
/// all used components need to follow this approach.
///
class IExtrapolationEngine
{
public:
  /// Virtual destructor
  virtual ~IExtrapolationEngine() {}

  /// Main extrapolation method, templated to chared/neutral
  ///
  /// @param ecCharged ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param bcheck is the boudnary check directive
  /// @todo shift to cell after splitting
  ///
  /// @return extrapolation code to indicate outcome
  virtual ExtrapolationCode
  extrapolate(ExCellCharged&       ecCharged,
              const Surface*       sf     = 0,
              const BoundaryCheck& bcheck = true) const = 0;

  /// Main extrapolation method, templated to chared/neutral
  ///
  /// @param ecNeutral ist he extrapolaiton cell
  /// @param sf is the (optional) destinaton surface
  /// @param bcheck is the boudnary check directive
  ///
  /// @return extrapolation code to indicate outcome
  virtual ExtrapolationCode
  extrapolate(ExCellNeutral&       ecNeutral,
              const Surface*       sf     = 0,
              const BoundaryCheck& bcheck = true) const = 0;

  /// define for which GeometrySignature this extrapolator is valid
  virtual GeometryType
  geometryType() const = 0;

protected:
  std::string m_sopPrefix;   /// prefix for screen output
  std::string m_sopPostfix;  /// prefix for screen output
};

}  // namespace