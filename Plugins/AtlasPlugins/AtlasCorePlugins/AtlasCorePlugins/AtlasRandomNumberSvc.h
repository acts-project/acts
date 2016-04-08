///////////////////////////////////////////////////////////////////
// AtlasRandomNumberSvc.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ATLASPLUGINS_RANDOMNUMBERSVC_H 
#define ATLASPLUGINS_RANDOMNUMBERSVC_H  1

// Gaudi include 
#include "GaudiKernel/ServiceHandle.h"

// Core module
#include "CoreInterfaces/ServiceBase.h"
#include "CoreInterfaces/IRandomNumberSvc.h"

// forward to Athena random number svc
#include "AthenaKernel/IAtRndmGenSvc.h"

// CLHEP include
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandGaussZiggurat.h"
#include "CLHEP/Random/RandLandau.h"
#include "CLHEP/Random/RandGamma.h"

namespace Acts {
  
  class AtlasRandomNumberSvc : public Acts::ServiceBase, virtual public IRandomNumberSvc {
        
  public:
    /** Constructor */
    AtlasRandomNumberSvc(const std::string& name, ISvcLocator* svc);
    
    /** Destructor */
    ~AtlasRandomNumberSvc();
    
    /** AlgTool initialize method */
    StatusCode initialize() final;

    /** AlgTool finalize method */
    StatusCode finalize() final;
    
    /** Query the interfaces **/
    StatusCode queryInterface( const InterfaceID& riid, void** ppvInterface );
    
    /** Implementation of the draw method */
    double draw(Acts::Distribution dist, double k=1.0, double lambda=1.0) const final {
      
      double rnd = 0.;
      
      switch (dist) {
	// the uniform case
	case Acts::Distribution::Flat : 
	  rnd = CLHEP::RandFlat::shoot(m_randomEngine);
	  break;
	// the gaussian case
	case Acts::Distribution::Gauss : 
	  rnd = CLHEP::RandGauss::shoot(m_randomEngine);
	  break;
	case Acts::Distribution::GaussZiggurat : 
	  rnd = CLHEP::RandGaussZiggurat::shoot(m_randomEngine);
	  break;
	// the Landau case
	case Acts::Distribution::Landau :
	  rnd = CLHEP::RandLandau::shoot(m_randomEngine);
	  break;
	// the Gamma case
	case Acts::Distribution::Gamma : 
	  rnd = CLHEP::RandGamma::shoot(m_randomEngine,k,lambda);
	  break; 
	default:
	  MSG_ERROR("Requested unknown distribution. Returning 0.");
      }
      return rnd;
    }
    
    
  private:
    /** Atlas Random Generator service  */
    ServiceHandle<IAtRndmGenSvc>                     m_rndGenSvc;
    
    /** Random engine  */
    CLHEP::HepRandomEngine*                          m_randomEngine;
       
    /** Name of the random number stream */
    std::string                                      m_randomEngineName;

    
  };
  
} // end of namespace

#endif // ATLASPLUGINS_RANDOMNUMBERSVC_H

