///////////////////////////////////////////////////////////////////
// IExtrapolationEngine.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_TRKEXINTERFACES_IEXTRAPOLATIONENGINE_H
#define ATS_TRKEXINTERFACES_IEXTRAPOLATIONENGINE_H

// Gaudi
#include "GaudiKernel/IService.h"
// Trk
#include "GeometryUtils/GeometrySignature.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"
 
namespace Ats {
  
  static const InterfaceID IID_IExtrapolationEngine("IExtrapolationEngine", 1, 0);

  typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
  typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

  /** @class IExtrapolationEngine

      Extrapolation engine interface for Charged and Neutral parameters,
      it serves as the Master extrapolation interface but also as the specialized
      extrapolation engine ones.

      The ExtrapolationEngine is desinged as a thread safe const-correct service,
      all used components need to follow this approach.

      @author Andreas Salzburger -at - cern.ch
  */

  class IExtrapolationEngine : virtual public IService {
     public:

       /** Virtual destructor */
       virtual ~IExtrapolationEngine(){}

       /** AlgTool interface methods */
       static const InterfaceID& interfaceID() { return IID_IExtrapolationEngine; }

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
      std::string                                         m_sopPrefix;            //!< prefix for screen output
      std::string                                         m_sopPostfix;           //!< prefix for screen output

  };


} // end of namespace

#endif // ATS_TRKEXINTERFACES_IEXTRAPOLATIONENGINE_H

