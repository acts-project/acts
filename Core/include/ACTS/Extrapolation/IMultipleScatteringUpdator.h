///////////////////////////////////////////////////////////////////
// IMultipleScatteringUpdator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_EXTRAPOLATIONINTERFACES_IMULTIPLESCATTERINGUPDATOR_H
#define ACTS_EXTRAPOLATIONINTERFACES_IMULTIPLESCATTERINGUPDATOR_H 1

#include "ACTS/EventData/ParticleDefinitions.h"

namespace Acts {

  class MaterialProperties;
   
  /**@class IMultipleScatteringUpdator
     Interface class IMultipleScatteringUpdator
     
     @author Andreas.Salzburger@cern.ch
    */
  class IMultipleScatteringUpdator {
      
     public:
       /**Virtual destructor*/
       virtual ~IMultipleScatteringUpdator(){}
       
       
      /** Calculate the sigma on theta introduced by multiple scatteringt */
      virtual double sigmaSquare(const MaterialProperties& mat,
                                 double p,
                                 double pathcorrection,
                                 ParticleType particle=pion,
                                 double deltaE=0.) const = 0;  
  };


} // end of namespace

#endif // ACTS_EXTRAPOLATIONINTERFACES_IMULTIPLESCATTERINGUPDATOR_H
