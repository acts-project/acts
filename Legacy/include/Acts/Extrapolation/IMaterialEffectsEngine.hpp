// This file is part of the Acts project.
//
// Copyright (C) 2016-2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IMaterialEffectsEngine.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"
#include "Acts/Extrapolation/MaterialUpdateMode.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

/// @class IMaterialEffectsEngine
///
/// Material effects engine interface for charged and neutral
/// (fast track simulation) ,
/// the update is alwyas on the:
///  - eCell.leadParmaeters && eCell.leadLayer
///  - if eCell.leadPameters == eCell.startParamters clone to new parameters
///    else : update the new parameters
class IMaterialEffectsEngine
{
public:
  /// Virtual destructor
  virtual ~IMaterialEffectsEngine() {}

  /// Public charged material effects interface
  ///
  /// @param ecCharged is the charged extrapolaiton cell
  /// @param msurface is the material surface (for curvilinear parameters)
  /// @param dir is the additional direction prescription
  /// @param matupstage is the update stage (pre/full/post)
  ///
  /// @return extrapolation code to indicate the progress
  virtual ExtrapolationCode
  handleMaterial(ExCellCharged&      ecCharged,
                 const Surface*      msurface   = nullptr,
                 PropDirection       dir        = alongMomentum,
                 MaterialUpdateStage matupstage = fullUpdate) const = 0;

  /// Public neutral material effects interface
  ///
  /// @param ecNeutral is the neutral extrapolaiton cell
  /// @param msurface is the material surface (for curvilinear parameters)
  /// @param dir is the additional direction prescription
  /// @param matupstage is the update stage (pre/full/post)
  ///
  /// @return extrapolation code to indicate the progress
  virtual ExtrapolationCode
  handleMaterial(ExCellNeutral&      ecNeutral,
                 const Surface*      msurface   = nullptr,
                 PropDirection       dir        = alongMomentum,
                 MaterialUpdateStage matupstage = fullUpdate) const = 0;

protected:
  std::string m_sopPrefix;   ///< prefix for screen output
  std::string m_sopPostfix;  ///< prefix for screen output
};

}  // end of namespace