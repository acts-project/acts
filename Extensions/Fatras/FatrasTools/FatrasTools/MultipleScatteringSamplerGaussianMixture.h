///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGaussianMixture.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H 1
 
// Gaudi
#include "GaudiKernel/ServiceHandle.h"
// Core module
#include "CoreInterfaces/AlgToolBase.h"
#include "CoreInterfaces/IRandomNumberSvc.h"
#include "Algebra/AlgebraDefinitions.h"
// EventData module
#include "EventDataUtils/ParticleHypothesis.h"
// Geometry module
#include "Material/MaterialProperties.h"
// Fatras module
#include "FatrasInterfaces/IMultipleScatteringSampler.h"

namespace Acts {
 
  /** @class MultipleScatteringSamplerGaussianMixture
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
 
  class MultipleScatteringSamplerGaussianMixture : public AlgToolBase, virtual public IMultipleScatteringSampler {
     
  public:
      /** AlgTool like constructor */
      MultipleScatteringSamplerGaussianMixture(const std::string&,const std::string&,const IInterface*);
     
      /**Virtual destructor*/
      virtual ~MultipleScatteringSamplerGaussianMixture();
     
      /** AlgTool initailize method.*/
      StatusCode initialize() final;
     
      /** AlgTool finalize method */
      StatusCode finalize() final;
     
      /** Calculate the theta introduced by multiple scattering,
       *          according to the RutherFord-Scott Formula           
       */
      double simTheta(const Acts::MaterialProperties& mat,
                      double p,
                      double pathcorrection,
                      Acts::ParticleHypothesis particle=Acts::pion) const final;
     
  private:
                 
      /** Random Generator service  */
      ServiceHandle<Acts::IRandomNumberSvc>      m_rndGenSvc;  
      
      /** boolean switch to include log term  */
      bool                                       m_log_include;   
      
      /** modifies the Fruehwirth/Regler model to fit with G4 */
      bool                                       m_optGaussianMixtureG4;  
      
      /** struct of Particle Masses */
      static Acts::ParticleMasses s_particleMasses;
     
      /** main factor of Rutherford-Scott formula */
      static double               s_main_RutherfordScott;  
      /** log factor of Rutherford-Scott formula */
      static double               s_log_RutherfordScott;   
                                  
      /** main factor for Rossi-Greisen formula */
      static double               s_main_RossiGreisen;     
      /** main factor for Rossi-Greisen formula */
      static double               s_log_RossiGreisen;      
      
      /** ========= Gaussian mixture model Fruehwirth, Regler Nucl. Inst. Methods A 456(2001) ========= */
      /** Gaussian mixture model: Sigma parameter a0 */
      static double         s_gausMixSigma1_a0;     
      /** Gaussian mixture model: Sigma parameter a1 */
      static double         s_gausMixSigma1_a1;     
      /** Gaussian mixture model: Sigma parameter a2 */
      static double         s_gausMixSigma1_a2;     
      
      /** Gaussian mixture model: Epsilon parameter a0 */
      static double         s_gausMixEpsilon_a0;     
      /** Gaussian mixture model: Epsilon parameter a1 */
      static double         s_gausMixEpsilon_a1;     
      /** Gaussian mixture model: Epsilon parameter a2 */
      static double         s_gausMixEpsilon_a2;     
     
      /** Gaussian mixture model: Epsilon parameter b0 */
      static double         s_gausMixEpsilon_b0;     
      /** Gaussian mixture model: Epsilon parameter b1 */
      static double         s_gausMixEpsilon_b1;     
      /** Gaussian mixture model: Epsilon parameter b2 */
      static double         s_gausMixEpsilon_b2;    
       
      /** projection factor to scale the projected angle out of the plane */
      static double               s_projectionFactor;      

     
    };
 
 
} // end of namespace


#endif // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGAUSSIANMIXTURE_H