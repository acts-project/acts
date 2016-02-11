///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H
#define ATS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H 1

// Gaudi
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/GaudiException.h"
// Core moudle
#include "CoreInterfaces/ServiceBase.h"
// Extrapolation module
#include "ExtrapolationInterfaces/INavigationEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "ExtrapolationInterfaces/IPropagationEngine.h"
#include "ExtrapolationInterfaces/IMaterialEffectsEngine.h"
// Geometry module
#include "GeometryInterfaces/ITrackingGeometrySvc.h"
#include "Volumes/BoundarySurface.h"
// EventData module
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"

namespace Ats {

  class TrackingGeometry;

  /** @class StaticNavigationEntine

      The static navigation engine for finding the next volume,
      propagate to the boundary, can be shared with other engines that have a static frame.

      @author Andreas.Salzburger -at- cern.ch
    */
  class StaticNavigationEngine : public ServiceBase, virtual public INavigationEngine {

      public:
        /** Constructor */
        StaticNavigationEngine(const std::string& name, ISvcLocator* svc);

        /** Destructor */
        ~StaticNavigationEngine();

        /** AlgTool initialize method */
        StatusCode initialize() final;

        /** AlgTool finalize method */
        StatusCode finalize() final;

        /** avoid method shaddowing */
        using INavigationEngine::resolveBoundary;
        using INavigationEngine::resolvePosition;

        /** resolve the boundary situation - for charged particles */
        ExtrapolationCode resolveBoundary(Ats::ExCellCharged& eCell, PropDirection dir=alongMomentum) const final;                                                                                          

        /** resolve the boundary situation - for neutral particles */
        ExtrapolationCode resolveBoundary(Ats::ExCellNeutral& eCelll, PropDirection dir=alongMomentum) const final;

        /** resolve the boundary situation - for charged particles */
        ExtrapolationCode resolvePosition(Ats::ExCellCharged& eCell, PropDirection dir=alongMomentum, bool noLoop=false) const final;          

        /** resolve the boundary situation - for neutral particles */
        ExtrapolationCode resolvePosition(Ats::ExCellNeutral& eCelll, PropDirection dir=alongMomentum, bool noLoop=false) const final;

        /** acces to tracking geometry */
        const TrackingGeometry& trackingGeometry() const throw (GaudiException);

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


        //!< retrieve TrackingGeometry
        StatusCode  updateTrackingGeometry() const;

        ServiceHandle<IPropagationEngine>                    m_propagationEngine;        //!< the used propagation engine
        ServiceHandle<IMaterialEffectsEngine>                m_materialEffectsEngine;    //!< the material effects updated

        ServiceHandle<ITrackingGeometrySvc>                  m_trackingGeometrySvc;       //!< ToolHandle to the TrackingGeometrySvc
        mutable const TrackingGeometry*                      m_trackingGeometry;          //!< the tracking geometry owned by the navigator
        std::string                                          m_trackingGeometryName;      //!< Name of the TrackingGeometry as given in Detector Store

    };

inline const Ats::TrackingGeometry& StaticNavigationEngine::trackingGeometry() const throw (GaudiException) {
    if (!m_trackingGeometry && updateTrackingGeometry().isFailure()){
        EX_MSG_FATAL("", "updateGeo", "", "Could not load TrackingGeometry with name '" << m_trackingGeometryName << "'. Aborting." );
        throw GaudiException("StaticNavigationEngine", "Problem with TrackingGeometry loading.", StatusCode::FAILURE);
    }
    return (*m_trackingGeometry);
}

} // end of namespace

//!< define the templated function
#include "StaticNavigationEngine.icc"

#endif // ATS_EXTRAPOLATIONENGINE_STATICNAVIGATIONENGINE_H

