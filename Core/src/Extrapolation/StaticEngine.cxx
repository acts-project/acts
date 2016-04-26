///////////////////////////////////////////////////////////////////
// StaticEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ACTS/Extrapolation/StaticEngine.h"

// constructor
Acts::StaticEngine::StaticEngine(const Acts::StaticEngine::Config& seConfig) :
  m_config()
{
  setConfiguration(seConfig);
}

// destructor
Acts::StaticEngine::~StaticEngine()
{}

// configuration
void Acts::StaticEngine::setConfiguration(const Acts::StaticEngine::Config& seConfig) {
  // steering of the screen outoput (SOP)
  IExtrapolationEngine::m_sopPrefix  = seConfig.prefix;
  IExtrapolationEngine::m_sopPostfix = seConfig.postfix;
  // copy the configuration 
  m_config = seConfig;
}   

// charged extrapolation /
Acts::ExtrapolationCode Acts::StaticEngine::extrapolate(ExCellCharged& ecCharged,
                                                      const Surface* sf,
                                                      const BoundaryCheck& bcheck) const
{ return extrapolateT<TrackParameters>(ecCharged,sf,ecCharged.propDirection,bcheck); }


// neutral extrapolation 
Acts::ExtrapolationCode Acts::StaticEngine::extrapolate(ExCellNeutral& ecNeutral,
                                                      const Surface* sf,
                                                      const BoundaryCheck& bcheck) const
{ return extrapolateT<NeutralParameters>(ecNeutral,sf,ecNeutral.propDirection,bcheck); }                                                  
