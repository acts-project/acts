///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ExtrapolationEngine/StaticNavigationEngine.h"

DECLARE_COMPONENT(Acts::StaticNavigationEngine)

// constructor
Acts::StaticNavigationEngine::StaticNavigationEngine(const std::string& name, ISvcLocator* svc) :
  Acts::ServiceBase(name, svc),
  m_propagationEngine("", name),
  m_materialEffectsEngine("", name),
  m_trackingGeometrySvc("", name),
  m_trackingGeometry(nullptr),
  m_trackingGeometryName("TrackingGeometry")
{
    // Services needed
    declareProperty("TrackingGeometrySvc"    , m_trackingGeometrySvc);
    declareProperty("PropagationEngine"      , m_propagationEngine);
    declareProperty("MaterialEffectsEngine"  , m_materialEffectsEngine);
    // The TrackingGeometry
    declareProperty("TrackingGeometry"       , m_trackingGeometryName);
    // steering of the screen outoput (SOP)
    declareProperty("OutputPrefix"           , m_sopPrefix);
    declareProperty("OutputPostfix"          , m_sopPostfix);
}

// destructor
Acts::StaticNavigationEngine::~StaticNavigationEngine()
{}


/** Query the interfaces. */
StatusCode Acts::StaticNavigationEngine::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
  if ( IID_INavigationEngine == riid )
    *ppvInterface = (INavigationEngine*)this;
  else  {
    // Interface is not directly available: try out a base class
    return Service::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}  

// the interface method initialize
StatusCode Acts::StaticNavigationEngine::initialize()
{
    MSG_DEBUG("initialize()" );
    // retrieve the tracking geometry servcie - crucial, abort when it can not be retrieved
    RETRIEVE_FATAL(m_trackingGeometrySvc);
    m_trackingGeometryName = m_trackingGeometrySvc->trackingGeometryName();
    // retrieve the propagation servcie - crucial, abort when it can not be retrieved
    RETRIEVE_FATAL(m_propagationEngine);
    // retrieve the material effects service do nothing if empty
    RETRIEVE_NONEMPTY_FATAL(m_materialEffectsEngine);
    // initialize done successfully
    return StatusCode::SUCCESS;
}

// the interface method finalize
StatusCode Acts::StaticNavigationEngine::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}

/** charged situation */
Acts::ExtrapolationCode Acts::StaticNavigationEngine::resolveBoundary(Acts::ExCellCharged& ecCharged, PropDirection dir) const
{ return resolveBoundaryT<Acts::TrackParameters>(ecCharged,dir); }

/** charged situation */
Acts::ExtrapolationCode Acts::StaticNavigationEngine::resolveBoundary(Acts::ExCellNeutral& ecNeutral, PropDirection dir) const
{ return resolveBoundaryT<Acts::NeutralParameters>(ecNeutral,dir); }

/** charged  */
Acts::ExtrapolationCode Acts::StaticNavigationEngine::resolvePosition(Acts::ExCellCharged& ecCharged, PropDirection dir, bool noLoop) const
{ return resolvePositionT<Acts::TrackParameters>(ecCharged,dir,noLoop); }

/** neutral */
Acts::ExtrapolationCode Acts::StaticNavigationEngine::resolvePosition(Acts::ExCellNeutral& ecNeutral, PropDirection dir, bool noLoop) const
{ return resolvePositionT<Acts::NeutralParameters>(ecNeutral,dir,noLoop); }

StatusCode Acts::StaticNavigationEngine::updateTrackingGeometry() const {
    // retrieve the TrackingGeometry from the particular TrackingGeometrySvc
    m_trackingGeometry = m_trackingGeometrySvc->trackingGeometry();
    if (!m_trackingGeometry) {
        EX_MSG_FATAL("", "updateGeom", "", "could not retrieve TrackingGeometry '" << m_trackingGeometryName  << "' from TrackingGeometryService." );
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}
