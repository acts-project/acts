// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialEffectsEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONENGINE_MATERIALEFFECTSENGINE_H
#define ACTS_EXTRAPOLATIONENGINE_MATERIALEFFECTSENGINE_H 1

#include "ACTS/Extrapolation/IMaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/MaterialUpdateMode.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Extrapolation/detail/MaterialInteraction.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

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
	std::shared_ptr<Logger>                               logger;
          bool          eLossCorrection;         //!< apply the energy loss correction
          bool          eLossMpv;                //!< apply the energy loss correction as most probable value
          bool          mscCorrection;           //!< apply the multiple (coulomb) scattering correction
          std::string   prefix;                  //!< screen output prefix
          std::string   postfix;                 //!< screen output postfix
          
          Config() :
	    logger(getDefaultLogger("MaterialEffectsEngine",Logging::INFO)),
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
      Config                 m_config;                //!< configuration struct
     
    private:
      const Logger& logger() const {return *m_config.logger;}
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
  inline MaterialEffectsEngine::Config MaterialEffectsEngine::getConfiguration() const { return m_config; }
  
} // end of namespace

#endif // ACTS_EXTRAPOLATIONENGINE_MATERIAKEFFECTSENGINE_H

