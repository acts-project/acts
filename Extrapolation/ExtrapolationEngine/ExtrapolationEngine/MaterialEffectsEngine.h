///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_TRKEXINTERFACES_MATERIAKEFFECTSENGINE_H
#define ATS_TRKEXINTERFACES_MATERIAKEFFECTSENGINE_H

// Ats
#include "CoreInterfaces/ServiceBase.h"
// Trk
#include "ExtrapolationInterfaces/IMaterialEffectsEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "ExtrapolationUtils/MaterialInteraction.h"
#include "ExtrapolationUtils/MaterialUpdateMode.h"
#include "EventDataUtils/PropDirection.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"
 
namespace Ats {
  
  class Layer;

  /** @class MaterialEffectsEngine

      Material effects engine interface for charged and neutral (fast track simulation) ,
      the update is alwyas on the:
       - eCell.leadParmaeters && eCell.leadLayer
       - if eCell.leadPameters == eCell.startParamters clone to new parameters
         else : update the new parameters

      @author Andreas Salzburger -at - cern.ch
  */
  class MaterialEffectsEngine : public Ats::ServiceBase, virtual public IMaterialEffectsEngine {
    public:

      /** Constructor */
      MaterialEffectsEngine(const std::string& name, ISvcLocator* svc);

      /** Destructor */
      ~MaterialEffectsEngine();

      /** AlgTool initialize method */
      StatusCode initialize() final;

      /** AlgTool finalize method */
      StatusCode finalize() final;

      /** charged extrapolation */
      ExtrapolationCode handleMaterial(ExCellCharged& ecCharged,
                                       PropDirection dir=alongMomentum,
                                       MaterialUpdateStage matupstage=fullUpdate) const final;

      /** neutral extrapolation - only for Fatras, dummy implementation here */
      ExtrapolationCode handleMaterial(ExCellNeutral& ecNeutral,
                                       PropDirection dir=alongMomentum,
                                       MaterialUpdateStage matupstage=fullUpdate) const final;

    protected:
      /** charged extrapolation - depending on the MaterialUpdateStage -
        it either manipulates or creates new parameters
        @TODO check for memory handling
        */
      const TrackParameters* updateTrackParameters(const Ats::TrackParameters& parameters,
                                                   Ats::ExCellCharged& eCell,
                                                   Ats::PropDirection dir,
                                                   Ats::MaterialUpdateStage matupstage) const; 
        
      MaterialInteraction                          m_interactionFormulae;     //!< the formulas concentrated
      ParticleMasses                               m_particleMasses;          //!< struct of Particle masses
      bool                                         m_eLossCorrection;         //!< apply the energy loss correction
      bool                                         m_eLossMpv;                //!< apply the energy loss correction as most probable value
      bool                                         m_mscCorrection;           //!< apply the multiple (coulomb) scattering correction

  };

} // end of namespace

#endif // ATS_TRKEXINTERFACES_MATERIAKEFFECTSENGINE_H

