///////////////////////////////////////////////////////////////////
// ExtrapolationEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Extrapolation/ExtrapolationEngine.h"

// constructor
Acts::ExtrapolationEngine::ExtrapolationEngine(const Acts::ExtrapolationEngine::Config& eeConfig) :
  m_eeConfig()
{
  setConfiguration(eeConfig);
}

// destructor
Acts::ExtrapolationEngine::~ExtrapolationEngine()
{}

// configuration
void Acts::ExtrapolationEngine::setConfiguration(const Acts::ExtrapolationEngine::Config& eeConfig) {
  // steering of the screen outoput (SOP)
  IExtrapolationEngine::m_sopPrefix  = eeConfig.prefix;
  IExtrapolationEngine::m_sopPostfix = eeConfig.postfix;
  // copy the configuration 
  m_eeConfig = eeConfig;
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
