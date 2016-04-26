///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H
#define ACTS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H 1

#include "ACTS/Extrapolation/INavigationEngine.h"
#include "ACTS/Extrapolation/ExtrapolationCell.h"
#include "ACTS/Extrapolation/IMaterialEffectsEngine.h"
#include "ACTS/Extrapolation/IPropagationEngine.h"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.h"
#include "ACTS/Volumes/BoundarySurface.h"
#include "ACTS/EventData/TrackParameters.h"
#include "ACTS/EventData/NeutralParameters.h"

namespace Acts {

  class TrackingGeometry;

  /** @class StaticNavigationEntine

      The static navigation engine for finding the next volume,
      propagate to the boundary, can be shared with other engines that have a static frame.

      @author Andreas.Salzburger -at- cern.ch
    */
  class StaticNavigationEngine : virtual public INavigationEngine {

      public:

        /** @struct Config 
            configuration struct for the StaticNavigationEngine
        */
        struct Config {
            
            std::shared_ptr<const IPropagationEngine>       propagationEngine;     //!< the used propagation engine
            std::shared_ptr<const IMaterialEffectsEngine>   materialEffectsEngine; //!< the material effects updated
            const TrackingGeometry*                         trackingGeometry;      //!< the tracking geometry used by the navigator
            
            std::string                                     prefix;                //!< output prefix
            std::string                                     postfix;               //!< output postfix
            
            Config() :
              trackingGeometry(nullptr),
              propagationEngine(nullptr),
              materialEffectsEngine(nullptr),
              prefix("[SN] - "),
              postfix(" - ")
            {}             
            
        };

        /** Constructor */
        StaticNavigationEngine(const Config& snConfig);

        /** Destructor */
        ~StaticNavigationEngine();

        /** avoid method shaddowing */
        using INavigationEngine::resolveBoundary;
        using INavigationEngine::resolvePosition;

        /** resolve the boundary situation - for charged particles */
        ExtrapolationCode resolveBoundary(ExCellCharged& eCell, PropDirection dir=alongMomentum) const final;                                                                                          

        /** resolve the boundary situation - for neutral particles */
        ExtrapolationCode resolveBoundary(ExCellNeutral& eCelll, PropDirection dir=alongMomentum) const final;

        /** resolve the boundary situation - for charged particles */
        ExtrapolationCode resolvePosition(ExCellCharged& eCell, PropDirection dir=alongMomentum, bool noLoop=false) const final;          

        /** resolve the boundary situation - for neutral particles */
        ExtrapolationCode resolvePosition(ExCellNeutral& eCelll, PropDirection dir=alongMomentum, bool noLoop=false) const final;
        
        /** Set configuration method */
        void setConfiguration(const Config& meConfig);
      
        /** Get configuration method */
        Config getConfiguration() const;

    protected:
        /** the configuration member of the static navigation engine */                                                     
        Config  m_snConfig; 
        
    private:
        /** resolve the boundary situation */
        template <class T> ExtrapolationCode resolveBoundaryT(ExtrapolationCell<T>& eCell,
                                                              PropDirection dir=alongMomentum) const;

        /** resolve position */
        template <class T> ExtrapolationCode resolvePositionT(ExtrapolationCell<T>& eCell, 
                                                              PropDirection dir=alongMomentum,
                                                              bool noLoop=false) const;

        /** deal with the boundary Surface - called by resolveBoundary */
        template <class T> ExtrapolationCode handleBoundaryT(ExtrapolationCell<T>& eCell,
                                                             const BoundarySurface<TrackingVolume>& bSurfaceTV,
                                                             PropDirection dir=alongMomentum,
                                                             bool stepout=false) const;

    };
    
    /** Return the configuration object */    
    inline StaticNavigationEngine::Config StaticNavigationEngine::getConfiguration() const { return m_snConfig; }

} // end of namespace

//!< define the templated function
#include "ACTS/Extrapolation/detail/StaticNavigationEngine.icc"

#endif // ACTS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H

