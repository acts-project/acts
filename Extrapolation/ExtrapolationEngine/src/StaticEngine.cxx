///////////////////////////////////////////////////////////////////
// StaticEngine.cxx, ATS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Trk include
#include "ExtrapolationEngine/StaticEngine.h"

DECLARE_COMPONENT(Ats::StaticEngine)

// constructor
Ats::StaticEngine::StaticEngine(const std::string& name, ISvcLocator* svc) :
  Ats::ServiceBase(name, svc),
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
Ats::StaticEngine::~StaticEngine()
{}

// the interface method initialize
StatusCode Ats::StaticEngine::initialize()
{
    MSG_DEBUG("initialize()" );
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
StatusCode Ats::StaticEngine::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}

/** charged extrapolation */
Ats::ExtrapolationCode Ats::StaticEngine::extrapolate(ExCellCharged& ecCharged,
                                                      const Surface* sf,
                                                      const BoundaryCheck& bcheck) const
{ return extrapolateT<TrackParameters>(ecCharged,sf,ecCharged.propDirection,bcheck); }


/** neutral extrapolation */
Ats::ExtrapolationCode Ats::StaticEngine::extrapolate(ExCellNeutral& ecNeutral,
                                                      const Surface* sf,
                                                      const BoundaryCheck& bcheck) const
{ return extrapolateT<NeutralParameters>(ecNeutral,sf,ecNeutral.propDirection,bcheck); }
