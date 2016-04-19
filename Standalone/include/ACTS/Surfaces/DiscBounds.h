///////////////////////////////////////////////////////////////////
// DiscBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACESDISCBOUNDS_H
#define ACTS_SURFACESDISCBOUNDS_H

// Trk included
#include "ACTS/Surfaces/SurfaceBounds.h"

namespace Acts {

  /**
   @class DiscBounds
   
   common base class for all bounds that are in a r/phi frame
    - simply introduced to avoid wrong bound assigments to surfaces
    
   @author Andreas.Salzburger@cern.ch
   */
        
      
  class DiscBounds : public SurfaceBounds {

    public:

      /** Default Constructor */
      DiscBounds() : SurfaceBounds() {}
      
      /** Destructor */
      virtual ~DiscBounds(){}
      
      /** Virtual Constructor */
      virtual DiscBounds* clone() const = 0;
      
  };
                              
} // end of namespace

#endif // ACTS_SURFACESDISCBOUNDS_H
