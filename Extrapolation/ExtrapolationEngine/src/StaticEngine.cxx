///////////////////////////////////////////////////////////////////
// StaticEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ExtrapolationEngine/StaticEngine.h"

DECLARE_SERVICE_FACTORY(Acts::StaticEngine)

// constructor
Acts::StaticEngine::StaticEngine(const std::string& name, ISvcLocator* svc) :
  Acts::ServiceBase(name, svc),
  m_propagationEngine("", name),          
  m_navigationEngine("", name),    
  m_materialEffectsEngine("", name)
{
    // The Tools needed
    declareProperty("PropagationEngine"                     , m_propagationEngine);
    declareProperty("NavigationEngine"                      , m_navigationEngine);
    declareProperty("MaterialEffectsEngine"                 , m_materialEffectsEngine);
    // steering of the screen outoput (SOP)
    declareProperty("OutputPrefix"                          , m_sopPrefix);
    declareProperty("OutputPostfix"                         , m_sopPostfix);
}

// destructor
Acts::StaticEngine::~StaticEngine()
{}

/** Query the interfaces. */
StatusCode Acts::StaticEngine::queryInterface(const InterfaceID& riid, void** ppvInterface)
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
StatusCode Acts::StaticEngine::initialize()
{
    MSG_DEBUG("initialize()" );
    //Service needs to be initialized
    if (!ServiceBase::initialize()) return StatusCode::FAILURE;
    // retrieve the propagation servcie - crucial, abort when it can not be retrieved
    RETRIEVE_FATAL(m_propagationEngine);
    // retrieve the navigation service   - crucial, abort when it can not be retrieved
    RETRIEVE_FATAL(m_navigationEngine);
    // retrieve the material effects service do nothing if empty
    RETRIEVE_NONEMPTY_FATAL(m_materialEffectsEngine);
    // initialize done successfully 
    return StatusCode::SUCCESS;
}    

// the interface method finalize
StatusCode Acts::StaticEngine::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}

/** charged extrapolation */
Acts::ExtrapolationCode Acts::StaticEngine::extrapolate(ExCellCharged& ecCharged,
                                                      const Surface* sf,
                                                      const BoundaryCheck& bcheck) const
{ return extrapolateT<TrackParameters>(ecCharged,sf,ecCharged.propDirection,bcheck); }


/** neutral extrapolation */
Acts::ExtrapolationCode Acts::StaticEngine::extrapolate(ExCellNeutral& ecNeutral,
                                                      const Surface* sf,
                                                      const BoundaryCheck& bcheck) const
{ return extrapolateT<NeutralParameters>(ecNeutral,sf,ecNeutral.propDirection,bcheck); }
