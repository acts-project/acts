///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ACTS/Extrapolation/StaticNavigationEngine.h"

// constructor
Acts::StaticNavigationEngine::StaticNavigationEngine(const Acts::StaticNavigationEngine::Config& snConfig) :
  m_config()
{    
  setConfiguration(snConfig);
}

// destructor
Acts::StaticNavigationEngine::~StaticNavigationEngine()
{}

// configuration
void Acts::StaticNavigationEngine::setConfiguration(const Acts::StaticNavigationEngine::Config& snConfig) {
  // steering of the screen outoput (SOP)
  INavigationEngine::m_sopPrefix  = snConfig.prefix;
  INavigationEngine::m_sopPostfix = snConfig.postfix;
  // copy the configuration 
  m_config = snConfig;
}

// charged situation 
Acts::ExtrapolationCode Acts::StaticNavigationEngine::resolveBoundary(Acts::ExCellCharged& ecCharged, PropDirection dir) const
{ return resolveBoundaryT<Acts::TrackParameters>(ecCharged,dir); }

// neutral situation
Acts::ExtrapolationCode Acts::StaticNavigationEngine::resolveBoundary(Acts::ExCellNeutral& ecNeutral, PropDirection dir) const
{ return resolveBoundaryT<Acts::NeutralParameters>(ecNeutral,dir); }

// chared situation
Acts::ExtrapolationCode Acts::StaticNavigationEngine::resolvePosition(Acts::ExCellCharged& ecCharged, PropDirection dir, bool noLoop) const
{ return resolvePositionT<Acts::TrackParameters>(ecCharged,dir,noLoop); }

