// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IMaterialEffectsEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H
#define ACTS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H 1

#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/MaterialUpdateMode.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"

namespace Acts {
  
  typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
  typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

  /** @class IMaterialEffectsEngine

      Material effects engine interface for charged and neutral (fast track simulation) ,
      the update is alwyas on the:
       - eCell.leadParmaeters && eCell.leadLayer
       - if eCell.leadPameters == eCell.startParamters clone to new parameters
         else : update the new parameters

  */
  class IMaterialEffectsEngine {
     public:

       /** Virtual destructor */
       virtual ~IMaterialEffectsEngine(){}

       /** charged extrapolation */
       virtual ExtrapolationCode handleMaterial(ExCellCharged& ecCharged,
                                                PropDirection dir=alongMomentum,
                                                MaterialUpdateStage matupstage=fullUpdate ) const = 0;


       /** neutral extrapolation */
       virtual ExtrapolationCode handleMaterial(ExCellNeutral& ecNeutral,
                                                PropDirection dir=alongMomentum,
                                                MaterialUpdateStage matupstage=fullUpdate) const = 0;

    protected:
      std::string                               m_sopPrefix;            //!< prefix for screen output
      std::string                               m_sopPostfix;           //!< prefix for screen output

  };


} // end of namespace

#endif // ACTS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H
