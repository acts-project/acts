///////////////////////////////////////////////////////////////////
// AtlasRandomNumberSvc.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>

// AtlasCorePlugins include
#include "AtlasCorePlugins/AtlasRandomNumberSvc.h"

DECLARE_COMPONENT(Acts::AtlasRandomNumberSvc)

// constructor
Acts::AtlasRandomNumberSvc::AtlasRandomNumberSvc(const std::string& name, ISvcLocator* svc):
  Acts::ServiceBase(name, svc),
  m_rndGenSvc("AtDSFMTGenSvc", name),
  m_randomEngine(nullptr),
  m_randomEngineName("FatrasRnd")
{
  // Random Number Service 
  declareProperty("RandomNumberService"  , m_rndGenSvc        , "Random number generator");
  declareProperty("RandomStreamName"     , m_randomEngineName , "Name of the random number stream");
}

// destructor
Acts::AtlasRandomNumberSvc::~AtlasRandomNumberSvc()
{}

/** Query the interfaces. */
StatusCode Acts::AtlasRandomNumberSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
  if ( IID_IRandomNumberSvc == riid )
    *ppvInterface = (IRandomNumberSvc*)this;
  else  {
    // Interface is not directly available: try out a base class
    return Service::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}


// the interface method initialize
StatusCode Acts::AtlasRandomNumberSvc::initialize()
{
    MSG_DEBUG("initialize()" );
    
    // get the random generator service - crucial, abort when it can not be retrieved 
    RETRIEVE_FATAL(m_rndGenSvc);
    
    //Get own engine with own seeds:
    m_randomEngine = m_rndGenSvc->GetEngine(m_randomEngineName);
    if (!m_randomEngine) {
      MSG_FATAL( "Could not get random engine '" << m_randomEngineName << "'" );
      return StatusCode::FAILURE;
   }
   
   // everything is fine!
   return StatusCode::SUCCESS;
}

// the interface method finalize
StatusCode Acts::AtlasRandomNumberSvc::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}