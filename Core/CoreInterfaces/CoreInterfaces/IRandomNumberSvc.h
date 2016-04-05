///////////////////////////////////////////////////////////////////
// IMagneticFieldSvc.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_CORE_IRANDOMNUMBERSVC_H
#define ACTS_CORE_IRANDOMNUMBERSVC_H 1

#ifdef ACTS_CORE_IRANDOMNUMBERSVC_PLUGIN
#include ACTS_CORE_IRANDOMNUMBERSVC_PLUGIN
#else 

// Gaudi includes
#include "GaudiKernel/IService.h"

namespace Acts {
  
  /** Declaration of the interface ID ( interface id, major version, minor version) */
  static const InterfaceID IID_IRandomNumberSvc("IRandomNumberSvc", 1, 0);
  
  
  /** @class IRandomNumberSvc
   * 
   * Random number service for drawing number with a specified distribution
   * 
   *  @author Noemi Calace -at- cern.ch
   */
  
  /** @enum Distribution
   * Distribution, enum for drawing number with a specified distribution
   */

  enum Distribution : unsigned int {
      Flat          = 0,
      Gauss         = 1,
      GaussZiggurat = 2,
      Landau        = 3,
      Gamma         = 4
  };
  
  class IRandomNumberSvc: virtual public IService {
    
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////
  public:
    /** Virtual destructor */
    virtual ~IRandomNumberSvc() {}
      
    /** Retrieve interface ID */
    static const InterfaceID& interfaceID() { return IID_IRandomNumberSvc; }

    /** draw the random number 
     * with a specified distribution */
    virtual double draw(Acts::Distribution, double k=1.0, double lambda=1.0) const = 0;
      
  };
}

#endif // ACTS_CORE_IRANDOMNUMBERSVC_PLUGIN

#endif //> ! ACTS_CORE_IRANDOMNUMBERSVC_H
