#ifndef ATS_PARAMETERDEFINITIONS_H
#define ATS_PARAMETERDEFINITIONS_H 1

#ifdef ATS_PARAMETER_DEFINITIONS_PLUGIN

#include ATS_PARAMETER_DEFINITIONS_PLUGIN

#else

namespace Ats
{
  enum ParDef : unsigned int
  {
    eLOC_1   = 0, ///< first coordinate in local surface frame
    eLOC_2   = 1, ///< second coordinate in local surface frame
    eLOC_R = eLOC_1,
    eLOC_PHI = eLOC_2,
    eLOC_RPHI = eLOC_1,
    eLOC_Z = eLOC_2,
    eLOC_X = eLOC_1,
    eLOC_Y = eLOC_2,
    eLOC_D0 = eLOC_1,
    eLOC_Z0 = eLOC_2,
    ePHI    = 2, ///< phi direction of momentum in global frame
    eTHETA  = 3, ///< theta direction of momentum in global frame
    eQOP    = 4,  ///< charge/momentum for charged tracks, for neutral tracks it is 1/momentum
    NGlobalPars
  };
}  // end of namespace Ats

#endif // ATS_PARAMETER_DEFINITIONS_PLUGIN

#endif //  ATS_PARAMETERDEFINITIONS_H
