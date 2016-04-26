///////////////////////////////////////////////////////////////////
// IMaterialEffectsEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H
#define ACTS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H 1

#include "ACTS/Extrapolation/ExtrapolationCell.h"
#include "ACTS/Extrapolation/MaterialUpdateMode.h"
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/EventData/TrackParameters.h"
#include "ACTS/EventData/NeutralParameters.h"

namespace Acts {
  
  typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
  typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

  /** @class IMaterialEffectsEngine

      Material effects engine interface for charged and neutral (fast track simulation) ,
      the update is alwyas on the:
       - eCell.leadParmaeters && eCell.leadLayer
       - if eCell.leadPameters == eCell.startParamters clone to new parameters
         else : update the new parameters

      @author Andreas Salzburger -at - cern.ch
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
