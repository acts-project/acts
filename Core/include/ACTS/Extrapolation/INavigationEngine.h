///////////////////////////////////////////////////////////////////
// INavigationEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONINTERFACES_INAVIGATIONENGINE_H
#define ACTS_EXTRAPOLATIONINTERFACES_INAVIGATIONENGINE_H

// Extrapolation module
#include "ExtrapolationUtils/ExtrapolationCell.h"
// EventData module
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"
 
namespace Acts {

  class TrackingGeometry;

  typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
  typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

  /** @class INavigationEngine

      Extrapolation engine interface for Charged and Neutral parameters,
      it serves as the Master extrapolation interface but also as the specialized
      extrapolation engine ones.

      @author Andreas Salzburger -at - cern.ch
  */

  class INavigationEngine : virtual public IService {

    public:

      /** Virtual destructor */
      virtual ~INavigationEngine(){}

      /** resolve the boundary situation - for charged particles */
      virtual ExtrapolationCode resolveBoundary(ExCellCharged& ecCell, PropDirection dir=alongMomentum) const = 0;

      /** resolve the boundary situation - for neutral particles */
      virtual ExtrapolationCode resolveBoundary(ExCellNeutral& enCell, PropDirection dir=alongMomentum) const = 0;

      /** resolve the position - for charged particles */
      virtual ExtrapolationCode resolvePosition(ExCellCharged& ecCell, PropDirection dir=alongMomentum, bool noLoop=false) const = 0;

      /** resolve the position - for neutral particles */
      virtual ExtrapolationCode resolvePosition(ExCellNeutral& enCell, PropDirection dir=alongMomentum, bool noLoop=false) const = 0;

      /** acces to tracking geometry */
      virtual const TrackingGeometry& trackingGeometry() const = 0;

    protected:
      //!< SCREEN output formatting  (SOP) - unify amongst extrapolation engines
      std::string m_sopPrefix;            //!< prefix for screen output
      std::string m_sopPostfix;           //!< prefix for screen output

  };


} // end of namespace

#endif // ACTS_EXTRAPOLATIONINTERFACES_INAVIGATIONENGINE_H

