///////////////////////////////////////////////////////////////////
// PhotonConversionSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_PHOTONCONVERSIONSAMPLER_H
#define ACTS_FATRASTOOLS_PHOTONCONVERSIONSAMPLER_H 1

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
// Core module
#include "CoreInterfaces/AlgToolBase.h"
#include "CoreInterfaces/IRandomNumberSvc.h"
#include "Algebra/AlgebraDefinitions.h"
// EventData module
#include "EventDataUtils/ParticleHypothesis.h"
#include "EventDataUtils/PdgToParticleHypothesis.h"
// Fatras module
#include "FatrasInterfaces/IPhotonConversionSampler.h"
// STL
#include <algorithm>

namespace Acts {
  
  class InteractionVertex;
  
  /**
   * @class PhotonConversionSampler
   * 
   * The PhotonConversionSampler called by the FatrasMaterialEffecsEngine
   * 
   * @author Sarka Todorova <Sarka.Todorova@cern.ch>
   * @author Noemi Calace   <Noemi.Calace@cern.ch>
   * 
  */
    class PhotonConversionSampler : public AlgToolBase, virtual public IPhotonConversionSampler {

  public:
    
    /** Constructor */
    PhotonConversionSampler(const std::string&,const std::string&,const IInterface*);
        
    /** Destructor */    
    ~PhotonConversionSampler();
    
    /** AlgTool interface methods */
    static const InterfaceID& interfaceID() { return IID_IPhotonConversionSampler; }
    
    /** AlgTool initialize method */
    StatusCode initialize() final;
        
    /** AlgTool finalize method */
    StatusCode finalize() final;
    
    /** processing of the presampled conversion on layer 
     * This methods return the photon's children
     */
    std::vector<Acts::InteractionVertex> doConversion(double time,
						      const Vector3D& position , 
						      const Vector3D& momentum) const final;
						      
    protected:
      
      double childEnergyFraction(double gammaMom) const;
      
      Vector3D childDirection(const Vector3D& gammaMom, double childE) const;
      
      std::vector<Acts::InteractionVertex> getChilds( double time, const Acts::Vector3D& vertex,
                                                      const Acts::Vector3D& photonMomentum,
                                                      double childEnergy, 
                                                      const Acts::Vector3D& childDirection,
                                                      Acts::ParticleHypothesis childType) const;
      /** helper functions for the Phi1/phi2 */
      double phi1(double delta) const;
      
      /** helper functions for the Phi1/phi2 */
      double phi2(double delta) const;
   
      /** Random Generator service  */
      ServiceHandle<Acts::IRandomNumberSvc>      m_rndGenSvc;  
      
      int                                        m_processCode;
      
      /** The cut from which on the child products are followed */
      double                                     m_minChildEnergy;
      double                                     m_childEnergyScaleFactor;
      
      /** struct of Particle Masses */
      static ParticleMasses                      s_particleMasses;
      static PdgToParticleHypothesis             s_pdgToHypo;
 
      /** Inverse fine structure constant */
      static double                              s_alpha;
      static double                              s_oneOverThree;
      
  };
  
  inline double PhotonConversionSampler::phi1(double delta) const {
    if (delta <= 1.)
      return 20.867 - 3.242 * delta  + 0.625*delta*delta;
    else
      return 21.12 - 4.184*log(delta+0.952);
  }
  
  inline double PhotonConversionSampler::phi2(double delta) const {
    if (delta <= 1.)
      return 20.209 - 1.930 * delta  - 0.086*delta*delta;
    return 21.12 - 4.184*log(delta+0.952);
  }
}

#endif // ACTS_FATRASTOOLS_PHOTONCONVERSIONSAMPLER_H