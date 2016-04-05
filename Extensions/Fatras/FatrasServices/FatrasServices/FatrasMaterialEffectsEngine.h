///////////////////////////////////////////////////////////////////
// FatrasMaterialEffectsEngine.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H
#define ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H 1

// Core module
#include "CoreInterfaces/ServiceBase.h"
#include "CoreInterfaces/IRandomNumberSvc.h"
// Athena & Gaudi includes
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
// Extrapolation module
#include "ExtrapolationInterfaces/IMaterialEffectsEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "ExtrapolationUtils/MaterialInteraction.h"
#include "ExtrapolationUtils/MaterialUpdateMode.h"
// EventData module
#include "EventDataUtils/PropDirection.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"

// Fatras module
#include "FatrasInterfaces/IEnergyLossSampler.h"
#include "FatrasInterfaces/IMultipleScatteringSampler.h"
#include "FatrasInterfaces/IPhotonConversionSampler.h"
#include "FatrasInterfaces/IHadronicInteractionSampler.h"
 
namespace Acts {
  
  /** @class FatrasMaterialEffectsEngine
   * 
   * Material effects engine interface for charged and neutral (fast track simulation) ,
   * the update is alwyas on the:
   * - eCell.leadParmaeters && eCell.leadLayer
   * - if eCell.leadPameters == eCell.startParamters clone to new parameters
   * else : update the new parameters
   * 
   * @author Andreas Salzburger <Andreas.Salzburger@cern.ch>
   * @author Noemi Calace       <Noemi.Calace@cern.ch>
   * 
   */
  class FatrasMaterialEffectsEngine : public Acts::ServiceBase, virtual public IMaterialEffectsEngine {
    public:

      /** Constructor */
      FatrasMaterialEffectsEngine(const std::string& name, ISvcLocator* svc);

      /** Destructor */
      ~FatrasMaterialEffectsEngine();

      /** AlgTool initialize method */
      StatusCode initialize() final;

      /** AlgTool finalize method */
      StatusCode finalize() final;

      /** Query the interfaces **/
      StatusCode queryInterface( const InterfaceID& riid, void** ppvInterface );

      /** charged extrapolation */
      ExtrapolationCode handleMaterial(ExCellCharged& ecCharged,
                                       PropDirection dir=alongMomentum,
                                       MaterialUpdateStage matupstage=fullUpdate) const final;

      /** neutral extrapolation - only for Fatras, dummy implementation here */
      ExtrapolationCode handleMaterial(ExCellNeutral& ecNeutral,
                                       PropDirection dir=alongMomentum,
                                       MaterialUpdateStage matupstage=fullUpdate) const final;

    protected:
      
      template <class T> ExtrapolationCode handleMaterialT( Acts::ExtrapolationCell<T>& eCell,
							    PropDirection dir=alongMomentum,
							    MaterialUpdateStage matupstage=fullUpdate) const;
							    
      ExtrapolationCode processMaterialOnLayer(ExCellCharged& ecCharged,
					       PropDirection dir,
					       float& mFraction) const;
      
      ExtrapolationCode processMaterialOnLayer(ExCellNeutral& ecNeutral,
					       PropDirection dir,
					       float& mFraction) const;
					       
      template <class T> ExtrapolationCode processMaterialOnLayerT( ExtrapolationCell<T>& eCell,
								    PropDirection dir,
								    float& mFraction) const;
								    
      const Acts::TrackParameters* updateTrackParameters(const Acts::TrackParameters& parameters,
							 Acts::ExCellCharged& eCell,
							 Acts::PropDirection dir,
							 double dX0,
							 double pathCorrection,
							 double mFraction) const;
							 
      const Acts::NeutralParameters* updateTrackParameters(const Acts::NeutralParameters& parameters,
							   Acts::ExCellNeutral& eCell,
							   Acts::PropDirection dir,
							   double dX0,
							   double pathCorrection,
							   double mFraction) const;
							   
      std::vector<Acts::InteractionVertex> interact(Acts::ExCellCharged& eCell, 
						    const Acts::Material&) const;

      std::vector<Acts::InteractionVertex> interact(Acts::ExCellNeutral& eCell, 
						    const Acts::Material&) const;
      
      void multipleScatteringUpdate(const Acts::TrackParameters& pars,
                                    ActsVectorD<5>& parameters,
                                    double simTheta, 
                                    double num_deltaPhi) const;
      
      void radiate( ActsVectorD<5> & parm ,
		    Acts::ExCellCharged& eCell, 
		    float pathLim, 
		    float mFr,
		    float refX) const;

      void collectBremPhoton(Acts::ExCellCharged& eCell,
			     double pElectron,
			     double gammaE,
			     const Vector3D& vertex,
			     Vector3D& particleDir) const;
      

      /** Random Generator service  */
      ServiceHandle<Acts::IRandomNumberSvc>            m_rndGenSvc;
      
      /** struct of Particle Masses */
      Acts::ParticleMasses                             m_particleMasses;
      
      /** IEnergyLossSampler */
      bool                                             m_doEnergyLoss;
      ToolHandle<Acts::IEnergyLossSampler>             m_energyLossSampler;
      
      /** Boolean switch for use of a dedicated eloss sampler */
      bool                                             m_dedicatedElectronEnergyLoss;  
      /** Pointer to the energy loss sampler - electrons */
      ToolHandle<Acts::IEnergyLossSampler>             m_elEnergyLossSampler;   
      
      /** Minimum momentum cut */
      double                                           m_minimumMomentum;
      
      /** Create brem photons flag */
      bool                                             m_createBremPhoton;
      /** Minimum momentum for brem photons */
      double                                           m_minimumBremPhotonMomentum;
      /** use the relativistic hertz dipole for brem photon radiation */
      bool                                             m_uniformHertzDipoleAngle;
      /** Process code for brem photons */
      int                                              m_bremProcessCode;
                  
      /** IMultipleScatteringSampler */
      bool                                             m_doMultipleScattering;
      ToolHandle<Acts::IMultipleScatteringSampler>     m_multipleScatteringSampler;
      /** describe deflection parametric/do real deflection */
      bool                                             m_parametricScattering;
      
      /** IPhotonConversionSampler */
      bool                                             m_doConversion;
      ToolHandle<Acts::IPhotonConversionSampler>       m_conversionSampler;
      
      /** IHadronicInteractionSampler */
      bool                                             m_doHadronicInteraction;
      ToolHandle<Acts::IHadronicInteractionSampler>    m_hadronicInteractionSampler;
      
      /** IDecaySampler */
      bool                                             m_doDecay;
      
      /** Do positrin annihilation */
      bool                                             m_doPositronAnnihilation;
      
      /** useful for the angle calculation of the brem photon */ 
      double                                           m_oneOverThree;
           
      /** useful for the multiple scattering calculation */ 
      double                                           m_projectionFactor;

    
  };

} // end of namespace

//!< define the templated function   
#include "FatrasMaterialEffectsEngine.icc" 

#endif // ACTS_FATRASSERVICES_FATRASMATERIALEFFECTSENGINE_H
