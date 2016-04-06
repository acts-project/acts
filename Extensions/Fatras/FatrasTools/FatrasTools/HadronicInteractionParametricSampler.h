///////////////////////////////////////////////////////////////////
// HadronicInteractionParametricSampler.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H
#define ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H 1

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
// Core module
#include "CoreInterfaces/AlgToolBase.h"
#include "CoreInterfaces/IRandomNumberSvc.h"
#include "Algebra/AlgebraDefinitions.h"
// EventData module
#include "EventDataUtils/ParticleHypothesis.h"
#include "EventDataUtils/PdgToParticleHypothesis.h"
// Fatras
#include "FatrasInterfaces/IHadronicInteractionSampler.h"
  
namespace Acts {
  
  class InteractionVertex;

  /** @class HadronicInteractionParametricSampler
   * 
   * Parametric implementation of nuclear interactions to be used
   * in Fatras. The parameterisation is gathered from Geant4.
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Carsten Magass     <Carsten.Magass@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   *                 
   */
   
  class HadronicInteractionParametricSampler : public AlgToolBase, virtual public IHadronicInteractionSampler {
  
  public:
    
    /** Constructor */
    HadronicInteractionParametricSampler(const std::string&,const std::string&,const IInterface*);
    
    /** Destructor */
    virtual ~HadronicInteractionParametricSampler();
    
    /** AlgTool initailize method.*/
    StatusCode initialize() final;
    
    /** AlgTool finalize method */
    StatusCode finalize() final;

    /** processing of the presampled nuclear interactions on layer
     * This method returns the particle's children
     */
    std::vector<Acts::InteractionVertex> doHadronicInteraction(double time,
							       const Acts::Vector3D& position, 
							       const Acts::Vector3D& momentum,
							       Acts::ParticleHypothesis particle=Acts::pion) const final;
 
   private:
     /** collect secondaries for layer material update */
     std::vector<Acts::InteractionVertex> getHadState( double time, double p,
						       const Acts::Vector3D& vertex,
						       const Acts::Vector3D& particleDir,
						       Acts::ParticleHypothesis particle ) const;
						       
     /** Random Generator service */
     ServiceHandle<Acts::IRandomNumberSvc>   m_rndGenSvc;
     
     /** process code */
     int                                     m_processCode;
     
     /** hadronic interaction setting */
     double                                  m_minimumHadOutEnergy;
     bool                                    m_cutChain;
              
     /** struct of Particle Masses */
     static ParticleMasses                   s_particleMasses;
     static PdgToParticleHypothesis          s_pdgToHypo;
           
   };
}

#endif // ACTS_FATRASTOOLS_HADRONICINTERACTIONPARAMETRICSAMPLER_H