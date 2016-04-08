///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerHighland.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H 1
 
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
// Extrapolation module
#include "ExtrapolationUtils/MaterialInteraction.h"

namespace Acts {
  
  /** @class MultipleScatteringSamplerHighland
   * 
   * The Formula used is the highland formula for the projected scattering angle :
   * 
   * @f$ \theta_{ms} = \frac{13.6MeV}{p}\cdot\sqrt{t/X_{0}}[1 + 0.038\ln(t/X_{0})] @f$
   * 
   * What is returned is the square of the expectation value of the deflection
   * @f$ < (\theta_ms)^2 > = \sigma_ms^2 @f$
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
 
  class MultipleScatteringSamplerHighland : public AlgToolBase, virtual public IMultipleScatteringSampler {
     
  public:
      /** AlgTool like constructor */
      MultipleScatteringSamplerHighland(const std::string&,const std::string&,const IInterface*);
     
      /**Virtual destructor*/
      virtual ~MultipleScatteringSamplerHighland();
     
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
      
      /** the formulas for multiple scattering evaluation */
      MaterialInteraction                        m_interactionFormulae;     
      
      /** boolean switch to include log term  */
      bool                                       m_log_include;           
      
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
                                  
      /** projection factor to scale the projected angle out of the plane */
      static double               s_projectionFactor;      

     
    };
 
 
} // end of namespace


#endif // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERHIGHLAND_H