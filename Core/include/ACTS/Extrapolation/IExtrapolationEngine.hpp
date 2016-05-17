///////////////////////////////////////////////////////////////////
// IExtrapolationEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONINTERFACES_IEXTRAPOLATIONENGINE_H
#define ACTS_EXTRAPOLATIONINTERFACES_IEXTRAPOLATIONENGINE_H 1

#include "ACTS/Utilities/GeometrySignature.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
 
namespace Acts {
  
  typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
  typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;
  
  class Surface;
  class BoundaryCheck;

  /** @class IExtrapolationEngine

      Extrapolation engine interface for Charged and Neutral parameters,
      it serves as the Master extrapolation interface but also as the specialized
      extrapolation engine ones.

      The ExtrapolationEngine is desinged as a thread safe const-correct service,
      all used components need to follow this approach.

      @author Andreas Salzburger -at - cern.ch
  */

  class IExtrapolationEngine {
     public:

       /** Virtual destructor */
       virtual ~IExtrapolationEngine(){}

       /** charged extrapolation */
       virtual ExtrapolationCode extrapolate(ExCellCharged& ecCharged,
                                             const Surface* sf = 0,
                                             const BoundaryCheck& bcheck = true) const = 0;


       /** neutral extrapolation */
       virtual ExtrapolationCode extrapolate(ExCellNeutral& ecNeutral,
                                             const Surface* sf = 0,
                                             const BoundaryCheck& bcheck = true) const = 0;


      /** define for which GeometrySignature this extrapolator is valid */
      virtual GeometryType geometryType() const = 0;

    protected:
      //!< SCREEN output formatting  (SOP) - unify amongst extrapolation engines
      std::string     m_sopPrefix;            //!< prefix for screen output
      std::string     m_sopPostfix;           //!< prefix for screen output

  };


} // end of namespace

#endif // ACTS_TRKEXINTERFACES_IEXTRAPOLATIONENGINE_H
