// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONENGING_STATICENGINE_H
#define ACTS_EXTRAPOLATIONENGING_STATICENGINE_H 1

#ifndef ACTS_EXTRAPOLATIONENINGE_OUTPUTHELPER 
#define ACTS_EXTRAPOLATIONENINGE_OUTPUTHELPER 1
#define OH_CHECKFOUND(object) ( object ? "found" : "not found")
#endif

#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/Extrapolation/IPropagationEngine.hpp"
#include "ACTS/Extrapolation/IMaterialEffectsEngine.hpp"
#include "ACTS/Extrapolation/INavigationEngine.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

    
  /** @class StaticEngine
    
    Extrapolation engine for static layer & volume setup.
    
    This engine relies on the fact that every position in a static layer setup can be uniquely associated to a layer (NavigationLayer or real physical layer),
    and thus to a voume at navigation level. The extrapolation process within a fully static setup is then realised as a step from layer to layer within a volume,
    and from volume to volume at a higher level.

    The entire code is written as a template in either charged or neutral parameters.  
    
  
  */
  class StaticEngine : virtual public IExtrapolationEngine {

      public:
          
        /** @enum ResolveLayerType 
            - use for code readability
        */
        enum ResolveLayerType {
          StartLayer                = 0,
          NavigationLayer           = 1,
          PassThroughLayer          = 2,
          SubStructureLayer         = 3,
          DestinationLayer          = 4,
          StartAndDestinationLayer  = 6,
          UndefinedLayer            = 5
        };          
        
        /** @struct Config 
            Configuration struct of the StaticEngine,
            holds the helper engines that can be configured
          */
        struct Config {
	std::shared_ptr<Logger>                               logger;
            
            std::shared_ptr<IPropagationEngine>     propagationEngine;        //!< the used propagation engine
            std::shared_ptr<INavigationEngine>      navigationEngine;         //!< the navigation engine to resolve the boundary
            std::shared_ptr<IMaterialEffectsEngine> materialEffectsEngine;    //!< the material effects updated  
            
            std::string                             prefix;                //!< output prefix
            std::string                             postfix;               //!< output postfix
            
            Config() :
	      logger(getDefaultLogger("StaticEngine",Logging::INFO)),
	      propagationEngine(nullptr),
              navigationEngine(nullptr),
              materialEffectsEngine(nullptr),
              prefix("[SE] - "),
              postfix(" - ")
            {}     
            
        };
        
        /** Constructor */
        StaticEngine(const Config& seConfig);
        
        /** Destructor */
        ~StaticEngine();
        
        using IExtrapolationEngine::extrapolate;
        
        /** charged extrapolation - public interface */
        ExtrapolationCode extrapolate(ExCellCharged& ecCharged,
                                      const Surface* sf = 0,
                                      const BoundaryCheck& bcheck = true) const final;


        /** neutral extrapolation - public interface */
        ExtrapolationCode extrapolate(ExCellNeutral& ecNeutral,
                                      const Surface* sf = 0,
                                      const BoundaryCheck& bcheck = true) const final;
                         
        /** define for which GeometrySignature this extrapolator is valid - this is GLOBAL */
        GeometryType geometryType() const final;                             

        /** Set configuration method */
        void setConfiguration(const Config& meConfig);
      
        /** Get configuration method */
        Config getConfiguration() const;  
        
    protected:
        /** Configuration struct */
        Config m_config;        
                         
     private:
      const Logger& logger() const {return *m_config.logger;}

    /** main loop extrapolation method */
        template <class T> ExtrapolationCode extrapolateT(ExtrapolationCell<T>& eCell,
                                                          const Surface* sf = 0,
                                                          PropDirection dir=alongMomentum,
                                                          const BoundaryCheck& bcheck = true) const;
                                                                  
        /** init Navigation for static setup */
        template <class T> ExtrapolationCode initNavigationT(ExtrapolationCell<T>& eCell,
                                                             const Surface* sf = 0,
                                                             PropDirection dir=alongMomentum,
                                                             const BoundaryCheck& bcheck = true) const;
                                                          
        /** main static layer handling */                                                  
        template <class T> ExtrapolationCode handleLayerT(ExtrapolationCell<T>& eCell,
                                                          const Surface* sf = 0,
                                                          PropDirection dir=alongMomentum,
                                                          const BoundaryCheck& bcheck = true) const;  

        /** main sub structure layer handling */                                                  
        template <class T> ExtrapolationCode resolveLayerT(ExtrapolationCell<T>& eCell,
                                                           const Acts::Surface* sf,
                                                           PropDirection dir=alongMomentum,
                                                           const BoundaryCheck& bcheck = true,
                                                           bool hasSubStructure = false,
                                                           bool isStartLayer = false,
                                                           bool isDestinationLayer =false) const;  
        /** handle the failure - as configured */
        template <class T> ExtrapolationCode handleReturnT(ExtrapolationCode eCode,
                                                           ExtrapolationCell<T>& eCell,
                                                           const Surface* sf = 0,
                                                           PropDirection dir=alongMomentum,
                                                           const BoundaryCheck& bcheck = true) const;
                                                           
            
    };

  /** Return the geomtry type */
  inline GeometryType StaticEngine::geometryType() const 
      { return Acts::Static; }
  
  /** Return the configuration object */    
  inline StaticEngine::Config StaticEngine::getConfiguration() const { return m_config; }
  


} // end of namespace

//!< define the templated function    
#include "ACTS/Extrapolation/detail/StaticEngine.icc"  

#endif // ACTS_EXTRAPOLATIONENGING_STATICENGINE_H

