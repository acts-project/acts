///////////////////////////////////////////////////////////////////
// StaticEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_EXTRAPOLATIONENGING_STATICENGINE_H
#define ATS_EXTRAPOLATIONENGING_STATICENGINE_H 1

#ifndef ATS_EXTRAPOLATIONENINGE_OUTPUTHELPER 
#define ATS_EXTRAPOLATIONENINGE_OUTPUTHELPER 1
#define OH_CHECKFOUND(object) ( object ? "found" : "not found")
#endif

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
// Core module
#include "CoreInterfaces/ServiceBase.h"
// Extrapolation module
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationInterfaces/IExtrapolationEngine.h"
#include "ExtrapolationInterfaces/IPropagationEngine.h"
#include "ExtrapolationInterfaces/IMaterialEffectsEngine.h"
#include "ExtrapolationInterfaces/INavigationEngine.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
// EventData module
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"

namespace Ats {

    
  /** @class StaticEngine
    
    Extrapolation engine for static layer & volume setup.
    
    This engine relies on the fact that every position in a static layer setup can be uniquely associated to a layer (NavigationLayer or real physical layer),
    and thus to a voume at navigation level. The extrapolation process within a fully static setup is then realised as a step from layer to layer within a volume,
    and from volume to volume at a higher level.

    The entire code is written as a template in either charged or neutral parameters.  
    
    @author Andreas.Salzburger -at- cern.ch
  
  */
  class StaticEngine : public Ats::ServiceBase, virtual public IExtrapolationEngine {

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
        
        /** Constructor */
        StaticEngine(const std::string& name, ISvcLocator* svc);
        
        /** Destructor */
        ~StaticEngine();

        /** AlgTool initialize method */
        StatusCode initialize() final;
        
        /** AlgTool finalize method */
        StatusCode finalize() final;
        
        /** Query the interfaces **/
        StatusCode queryInterface( const InterfaceID& riid, void** ppvInterface );
        
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
                         
     private:
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
                                                           const Ats::Surface* sf,
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
                                                           
        ServiceHandle<IPropagationEngine>      m_propagationEngine;        //!< the used propagation engine
        ServiceHandle<INavigationEngine>       m_navigationEngine;         //!< the navigation engine to resolve the boundary
        ServiceHandle<IMaterialEffectsEngine>  m_materialEffectsEngine;    //!< the material effects updated
            
    };

  inline GeometryType StaticEngine::geometryType() const 
      { return Ats::Static; }


} // end of namespace

//!< define the templated function    
#include "StaticEngine.icc"  

#endif // ATS_EXTRAPOLATIONENGING_STATICENGINE_H

