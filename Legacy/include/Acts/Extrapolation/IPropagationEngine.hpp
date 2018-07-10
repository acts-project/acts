// This file is part of the Acts project.
//
// Copyright (C) 2016-2018  Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IPropagationEngine.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"

namespace Acts {

typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

/// @class IPropagationEngine
///
/// A propagation engine wrapping the propagator algtool it respects the path
/// limit
/// to stop particles if needed.
///
/// If the propagation is successful to the surface it will return
/// SuccessfulDestination,
/// the parameters will be attached to the ExtrapolationCell as leadParameters,
/// such that the engine can chose.
///
/// It also wraps the MultiTrackParameters (@todo do this actually)
class IPropagationEngine
{
public:
  /// Virtual destructor
  virtual ~IPropagationEngine() {}
  /// Main Charged extrapolation method
  ///
  /// @param ecCell is the charged extrapolation cell
  /// @param sf is the destination surface
  /// @param dir is the additional direction prescription
  /// @param purpose sets the fill mode in to the ExtrapolationCache
  /// @param bcheck is the boundary check prescription
  /// @param returnCurvilinear is a boolean switch to not collapse onto the
  ///        surface frame but stay in curviliear coordinates
  ///
  /// @return possible return codes :
  ///  - SuccessPathLimit (path limit reached)
  ///  - SucessDestination (surface hit, only when finalPropagation == true)
  ///  - InProgress (surface hit, when finalPropagation == false)
  ///  - Recovered (surface not hit, leadParameters stay untouched)
  virtual ExtrapolationCode
  propagate(ExCellCharged&                        ecCell,
            const Surface&                        sf,
            NavigationDirection                   dir = forward,
            std::vector<ExtrapolationMode::eMode> purpose
            = {ExtrapolationMode::Destination},
            const BoundaryCheck& bcheck            = true,
            bool                 returnCurvilinear = true) const = 0;

  /// Main Neutral extrapolation method
  ///
  /// @param enCell is the neutral extrapolation cell
  /// @param sf is the destination surface
  /// @param dir is the additional direction prescription
  /// @param purpose sets the fill mode in to the ExtrapolationCache
  /// @param bcheck is the boundary check prescription
  /// @param returnCurvilinear is a boolean switch to not collapse onto the
  ///        surface frame but stay in curviliear coordinates
  ///
  /// @return possible return codes :
  ///  - SuccessPathLimit (path limit reached)
  ///  - SucessDestination (surface hit, only when finalPropagation == true)
  ///  - InProgress (surface hit, when finalPropagation == false)
  ///  - Recovered (surface not hit, leadParameters stay untouched)
  virtual ExtrapolationCode
  propagate(ExCellNeutral&                        enCell,
            const Surface&                        sf,
            NavigationDirection                   dir = forward,
            std::vector<ExtrapolationMode::eMode> purpose
            = {ExtrapolationMode::Destination},
            const BoundaryCheck& bcheck            = true,
            bool                 returnCurvilinear = true) const = 0;

protected:
  std::string m_sopPrefix;   ///< prefix for screen output
  std::string m_sopPostfix;  ///< prefix for screen output
};

}  // namespace