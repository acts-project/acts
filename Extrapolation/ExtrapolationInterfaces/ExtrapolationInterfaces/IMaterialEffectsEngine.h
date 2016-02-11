///////////////////////////////////////////////////////////////////
// IMaterialEffectsEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H
#define ATS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H 1

// Gaudi
#include "GaudiKernel/IService.h"
// Extrapolation module
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "ExtrapolationUtils/MaterialUpdateMode.h"
// EventData module
#include "EventDataUtils/PropDirection.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"

namespace Ats {
  
  static const InterfaceID IID_IMaterialEffectsEngine("IMaterialEffectsEngine", 1, 0);

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
  class IMaterialEffectsEngine : virtual public IService {
     public:

       /** Virtual destructor */
       virtual ~IMaterialEffectsEngine(){}

       /** AlgTool interface methods */
       static const InterfaceID& interfaceID() { return IID_IMaterialEffectsEngine; }

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

#endif // ATS_EXTRAPOLATIONINTERFACES_IMATERIAKEFFECTSENGINE_H
