///////////////////////////////////////////////////////////////////
// MultipleScatteringSamplerGeneralMixture.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGENERALMIXTURE_H
#define ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGENERALMIXTURE_H 1
 
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
  
   /** @class MultipleScatteringSamplerGeneralMixture
   *
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
 
  class MultipleScatteringSamplerGeneralMixture : public AlgToolBase, virtual public IMultipleScatteringSampler {
     
  public:
      /** AlgTool like constructor */
      MultipleScatteringSamplerGeneralMixture(const std::string&,const std::string&,const IInterface*);
     
      /**Virtual destructor*/
      virtual ~MultipleScatteringSamplerGeneralMixture();
     
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
                
      /** struct of Particle Masses */
      static Acts::ParticleMasses s_particleMasses;
     
      /** main factor for Rossi-Greisen formula */
      static double               s_main_RossiGreisen;     
      /** main factor for Rossi-Greisen formula */
      static double               s_log_RossiGreisen;      
      
      /** ========= General mixture model Fruehwirth, M. Liendl. Comp. Phys. Comm. 141 (2001) 230-246 ========= */
      
      /** General mixture model: Scaling factor */
      static double         s_genMixScale;          
      
      /** General mixture model: get parameters for single gaussian simulation */
      double*    getGaussian(double beta, double p,double dOverX0, double scale) const;
      /** General mixture model: get parameters for gaussian mixture */
      double*    getGaussmix(double beta, double p,double dOverX0,double Z, double scale) const;
      /** General mixture model: get parameters for semi-gaussian mixture */
      double*    getSemigauss(double beta,double p,double dOverX0,double Z, double scale) const;
      /** General mixture model: simulate semi-gaussian mixture */
      double    simGaussmix(double* scattering_params) const;
      /** General mixture model: simulate gaussian mixture */
      double    simSemigauss(double* scattering_params) const;
       
      /** projection factor to scale the projected angle out of the plane */
      static double               s_projectionFactor;      

     
    };
 
 
} // end of namespace


#endif // ACTS_FATRASTOOLS_MULTIPLESCATTERINGSAMPLERGENERALMIXTURE_H