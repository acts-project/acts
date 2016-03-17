///////////////////////////////////////////////////////////////////
// ExtrapolationEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ExtrapolationEngine/ExtrapolationEngine.h"
// Geometry module
#include "GeometryUtils/GeometrySignature.h"


DECLARE_COMPONENT(Acts::ExtrapolationEngine)

// constructor
Acts::ExtrapolationEngine::ExtrapolationEngine(const std::string& name, ISvcLocator* svc)
: Acts::ServiceBase(name, svc),
  m_trackingGeometry(nullptr),
  m_trackingGeometrySvc("TrackingGeometrySvc", name),
  m_trackingGeometryName("TrackingGeometry"),
  m_extrapolationEngines(name),
  m_propagationEngine("", name),    
  m_navigationEngine("", name),
  m_forceSearchInit(false)
{
    // Geometry retrieval
    declareProperty("TrackingGeometrySvc"                   , m_trackingGeometrySvc);
    // Extrapolation Engine retrieval 
    declareProperty("ExtrapolationEngines"                  , m_extrapolationEngines);    
    // The Tools needed
    declareProperty("PropagationEngine"                     , m_propagationEngine);
    declareProperty("NavigationEngine"                      , m_navigationEngine);
    // steering of the screen outoput (SOP)
    declareProperty("OutputPrefix"                          , m_sopPrefix);
    declareProperty("OutputPostfix"                         , m_sopPostfix);
    // the properties to be given 
    declareProperty("ForceSearchAtInit"                     , m_forceSearchInit);
}

// destructor
Acts::ExtrapolationEngine::~ExtrapolationEngine()
{}

/** Query the interfaces. */
StatusCode Acts::ExtrapolationEngine::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
  if ( IID_IExtrapolationEngine == riid )
    *ppvInterface = (IExtrapolationEngine*)this;
  else  {
    // Interface is not directly available: try out a base class
    return Service::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}
    
// the interface method initialize
StatusCode Acts::ExtrapolationEngine::initialize()
{            
    MSG_DEBUG("initialize()");
    // retrieve the tracking geometry servcie - crucial, abort when it can not be retrieved
    RETRIEVE_FATAL(m_trackingGeometrySvc);
    m_trackingGeometryName = m_trackingGeometrySvc->trackingGeometryName();    
    // retriveve the extrapolation engines - crucial, abort when they can not be retrieved
    RETRIEVE_FATAL(m_extrapolationEngines);
    EX_MSG_DEBUG( "", "initialize", "", "Successfully retrieved " << m_extrapolationEngines.size() << " ExtrapolationEngines. Ordering them now." );
    m_eeAccessor = std::vector<const Acts::IExtrapolationEngine*>(int(Acts::NumberOfGeometryTypes), (const Acts::IExtrapolationEngine*)nullptr);
    for (auto& ee : m_extrapolationEngines){
        EX_MSG_DEBUG( "", "initialize", "", "Registering " << ee->name() << " - for GeometryType : "  << ee->geometryType() );
        m_eeAccessor[ee->geometryType()] = (&*ee);
    }
    // retrive a propagation engine for initialization - crucial, abort when they can not be retrieved
    RETRIEVE_FATAL(m_propagationEngine);
    // retrieve a navigation engine - crucial, abort when they can not be retrieved
    RETRIEVE_FATAL(m_navigationEngine);
    // everything was successfully retrieved
    return StatusCode::SUCCESS;
}    

// the interface method finalize
StatusCode Acts::ExtrapolationEngine::finalize()
{    
    MSG_DEBUG("finalize()");
    return StatusCode::SUCCESS;
}

/** charged extrapolation */
Acts::ExtrapolationCode Acts::ExtrapolationEngine::extrapolate(ExCellCharged& ecCharged,
                                                        const Surface* sf,
                                                        const BoundaryCheck& bcheck) const
{ return extrapolateT<TrackParameters>(ecCharged,sf,ecCharged.propDirection,bcheck); }


/** neutral extrapolation */
Acts::ExtrapolationCode Acts::ExtrapolationEngine::extrapolate(ExCellNeutral& ecNeutral,
                                                        const Surface* sf,
                                                        const BoundaryCheck& bcheck) const
{ return extrapolateT<NeutralParameters>(ecNeutral,sf,ecNeutral.propDirection,bcheck); }
                                           

StatusCode Acts::ExtrapolationEngine::updateTrackingGeometry() const {
    // retrieve the TrackingGeometry from the particular TrackingGeometrySvc
    m_trackingGeometry = m_trackingGeometrySvc->trackingGeometry();
    if (!m_trackingGeometry) {
        EX_MSG_FATAL("", "updateGeom", "", "could not retrieve TrackingGeometry '" << m_trackingGeometryName  << "' from TrackingGeometryService." );
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}


