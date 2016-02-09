///////////////////////////////////////////////////////////////////
// IMultipleScatteringUpdator.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef TRKEXINTERFACES_IMULTIPLESCATTERINGUPDATOR_H
#define TRKEXINTERFACES_IMULTIPLESCATTERINGUPDATOR_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// Trk
#include "EventDataUtils/ParticleHypothesis.h"

namespace Ats {

  class MaterialProperties;
  
  /** Interface ID for IMultipleScatteringUpdator*/  
  static const InterfaceID IID_IMultipleScatteringUpdator("IMultipleScatteringUpdator", 1, 0);
  
  /**@class IMultipleScatteringUpdator
     Interface class IMultipleScatteringUpdator
     
     @author Andreas.Salzburger@cern.ch
    */
  class IMultipleScatteringUpdator : virtual public IAlgTool {
      
     public:
       /**Virtual destructor*/
       virtual ~IMultipleScatteringUpdator(){}
       
       /** AlgTool and IAlgTool interface methods */
       static const InterfaceID& interfaceID() { return IID_IMultipleScatteringUpdator; };
       
      /** Calculate the sigma on theta introduced by multiple scatteringt  
      */
      virtual double sigmaSquare(const MaterialProperties& mat,
                                 double p,
                                 double pathcorrection,
                                 ParticleHypothesis particle=pion,
                                 double deltaE=0.) const = 0;

     /** Validation Action:
          Can be implemented optionally, outside access to internal validation steps */
       virtual void validationAction() const {}
  
  };


} // end of namespace


#endif // TRKEXINTERFACES_IMULTIPLESCATTERINGUPDATOR_H


