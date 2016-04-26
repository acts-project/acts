///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONENGINE_MATERIALEFFECTSENGINE_H
#define ACTS_EXTRAPOLATIONENGINE_MATERIALEFFECTSENGINE_H 1

#include "ACTS/Extrapolation/IMaterialEffectsEngine.h"
#include "ACTS/Extrapolation/ExtrapolationCell.h"
#include "ACTS/Extrapolation/MaterialUpdateMode.h"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.h"
#include "ACTS/Extrapolation/detail/MaterialInteraction.h"
#include "ACTS/EventData/TrackParameters.h"
#include "ACTS/EventData/NeutralParameters.h"
#include "ACTS/Utilities/Definitions.h"
 
namespace Acts {
  
  class Layer;

  /** @class MaterialEffectsEngine

      Material effects engine interface for charged and neutral (fast track simulation) ,
      the update is alwyas on the:
       - eCell.leadParmaeters && eCell.leadLayer
       - if eCell.leadPameters == eCell.startParamters clone to new parameters
         else : update the new parameters

      @author Andreas Salzburger -at - cern.ch
  */
  class MaterialEffectsEngine : virtual public IMaterialEffectsEngine {

    public:
      /** @struct Config 
          Configuration struct for the MaterialEffectsEngine
        */
      struct Config {
          
          bool          eLossCorrection;         //!< apply the energy loss correction
          bool          eLossMpv;                //!< apply the energy loss correction as most probable value
          bool          mscCorrection;           //!< apply the multiple (coulomb) scattering correction
          std::string   prefix;                  //!< screen output prefix
          std::string   postfix;                 //!< screen output postfix
          
          Config() :
            eLossCorrection(true),
            eLossMpv(true),        
            mscCorrection(true),
            prefix("[ME] - "),
            postfix(" - ") 
         {}         
          
      };        

      /** Constructor */
      MaterialEffectsEngine(const Config& meConfig);

      /** Destructor */
      ~MaterialEffectsEngine();

      /** charged extrapolation */
      ExtrapolationCode handleMaterial(ExCellCharged& ecCharged,
                                       PropDirection dir=alongMomentum,
                                       MaterialUpdateStage matupstage=fullUpdate) const final;

      /** neutral extrapolation - only for Fatras, dummy implementation here */
      ExtrapolationCode handleMaterial(ExCellNeutral& ecNeutral,
                                       PropDirection dir=alongMomentum,
                                       MaterialUpdateStage matupstage=fullUpdate) const final;
                                       
                                       
      /** Set configuration method */
      void setConfiguration(const Config& meConfig);
      
      /** Get configuration method */
      Config getConfiguration() const;                                      

    protected:
      Config                 m_meConfig;                //!< configuration struct
     
    private:
      /** charged extrapolation - depending on the MaterialUpdateStage -
        it either manipulates or creates new parameters
        @TODO check for memory handling
        */
      const TrackParameters* updateTrackParameters(const Acts::TrackParameters& parameters,
                                                   Acts::ExCellCharged& eCell,
                                                   Acts::PropDirection dir,
                                                   Acts::MaterialUpdateStage matupstage) const; 
        
      MaterialInteraction    m_interactionFormulae;     //!< the formulas concentrated
      ParticleMasses         m_particleMasses;          //!< struct of Particle masses


  };
  
  /** Return the configuration object */    
  inline MaterialEffectsEngine::Config MaterialEffectsEngine::getConfiguration() const { return m_meConfig; }
  
} // end of namespace

#endif // ACTS_EXTRAPOLATIONENGINE_MATERIAKEFFECTSENGINE_H

