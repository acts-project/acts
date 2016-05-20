/////////////////////////////////////////////////////////////////////////////////
//  Header file for class RungeKuttaEngine, ACTS project
/////////////////////////////////////////////////////////////////////////////////

#ifndef ACTS_RUNGEKUTTAENGINE_RUNGEKUTAENGINE_H
#define ACTS_RUNGEKUTTAENGINE_RUNGEKUTAENGINE_H 1

#include "ACTS/Extrapolation/IPropagationEngine.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/detail/ExtrapolationMacros.hpp"
#include "ACTS/Extrapolation/detail/RungeKuttaUtils.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/ConeSurface.hpp"
#include "ACTS/Surfaces/BoundaryCheck.hpp"
#include "ACTS/MagneticField/IMagneticFieldSvc.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {

  /**
    @struct PropagationCache 
     Helper struct to allow state-less propagation.
    */
    
  struct PropagationCache { 
      
      // configuration
      double                    direction;
      BoundaryCheck             boundaryCheck;
      bool                      returnCurvilinear;
      bool                      useJacobian;
      double                    step;
      double                    maxPathLength;
      bool                      maxPathLimit;
      bool                      mcondition;
      bool                      needgradient;  
      bool                      newfield;
      // internal parameters to be used
      double                    field[3] = {0.,0.,0.};
      double                    pVector[64];
      // result
      double                    parameters[5] = {0.,0.,0.,0.,0.};
      ActsSymMatrixD<5>*        covariance; 
      double                    jacobian[25];
      
      
      PropagationCache() :
       direction(alongMomentum),
       boundaryCheck(true),
       returnCurvilinear(false),
       useJacobian(false),
       step(0.),
       maxPathLength(0.),
       maxPathLimit(false),
       mcondition(false),
       needgradient(false), 
       newfield(true),
       covariance(nullptr)
      {}
  
  };

  class Surface;

  /**
  @class RungeKuttaEngine
    
    RungeKuttaEngine is algorithm for track parameters propagation through
    magnetic field with or without jacobian of transformation. This algorithm
    contains three steps.
    
    1.The first step of the algorithm is track parameters transformation from
      local presentation for given start surface to global Runge Kutta coordinates.
    
    2.The second step is propagation through magnetic field with or without jacobian.
    
    3.Third step is transformation from global Runge Kutta presentation to local
      presentation of given output surface.
     
    
      AtaPlane    AtaStraightLine      AtaDisc       AtaCylinder      Perigee
         |               |               |               |              |
         |               |               |               |              |
         V               V               V               V              V 
         ----------------------------------------------------------------- 
                                         |              Local->Global transformation
                                         V
                      Global position (Runge Kutta presentation)
                                         |
                                         |
                   Propagation to next surface with or without jacobian
                             using Nystroem algorithm 
                 (See Handbook Net. Bur. of Standards, procedure 25.5.20)
                                         |
                                         V              Global->Local transformation
         ----------------------------------------------------------------
         |               |               |               |              |
         |               |               |               |              |
         V               V               V               V              V
     PlaneSurface StraightLineSurface DiscSurface CylinderSurface PerigeeSurface 
    
    For propagation using Runge Kutta method we use global coordinate, direction,
    inverse momentum and Jacobian of transformation. All this parameters we save 
    in array P[42] called pVector
  
                     /dL0    /dL1    /dPhi   /dThe   /dCM
    X  ->P[0]  dX /   P[ 7]   P[14]   P[21]   P[28]   P[35]  
    Y  ->P[1]  dY /   P[ 8]   P[15]   P[22]   P[29]   P[36]  
    Z  ->P[2]  dZ /   P[ 9]   P[16]   P[23]   P[30]   P[37]   
    Ax ->P[3]  dAx/   P[10]   P[17]   P[24]   P[31]   P[38]  
    Ay ->P[4]  dAy/   P[11]   P[18]   P[25]   P[32]   P[39]  
    Az ->P[5]  dAz/   P[12]   P[19]   P[26]   P[33]   P[40]  
    CM ->P[6]  dCM/   P[13]   P[20]   P[27]   P[34]   P[41] 
    
    where 
         in case local presentation 
    
         L0  - first  local coordinate  (surface dependent)
         L1  - second local coordinate  (surface dependent)
         Phi - Azimuthal angle
         The - Polar     angle
         CM  - charge/momentum
    
         in case global presentation
    
         X   - global x-coordinate        = surface dependent
         Y   - global y-coordinate        = surface dependent
         Z   - global z-coordinate        = sutface dependent
         Ax  - direction cosine to x-axis = Sin(The)*Cos(Phi)
         Ay  - direction cosine to y-axis = Sin(The)*Sin(Phi)
         Az  - direction cosine to z-axis = Cos(The)
         CM  - charge/momentum            = local CM
    
    Comment: 
         if pointer to const *  = 0 algorithm will propagate track 
         parameters and jacobian of transformation according straight line model
    
    
    @author Igor.Gavrilenko@cern.ch (adapted to ATS by Andreas.Salzburger -at- cern.ch)   
  */

  class RungeKuttaEngine : virtual public IPropagationEngine {
      
    public:
        
      /** @struct Config
          Configuration struct for the RungeKuttaEngine
        
          @TODO explain parametr meanings (input from Igor needed)
        */
      struct Config {
	std::shared_ptr<Logger>                               logger;
    
        std::shared_ptr<IMagneticFieldSvc>   fieldService  ;  //!< the field service 
        double                               dlt           ;  //!< accuracy parameter
        double                               helixStep     ;  //!< max step whith helix model
        double                               straightStep  ;  //!< max step whith srtaight line model
        double                               maxPathLength ;  //!< max overal path length 
        bool                                 usegradient   ;  //!< use magnetif field gradient
        std::string                          prefix        ;  //!< screen output prefix
        std::string                          postfix       ;  //!< screen output postfix
    
        Config() :
	  logger(getDefaultLogger("RungeKuttaEngine",Logging::INFO)),
          fieldService(nullptr),
          dlt(0.000200),           
          helixStep(1.),     
          straightStep(0.01),  
          maxPathLength(25000.),
          usegradient(false),
          prefix("[RK] - "),
          postfix(" - ")
         {}
           
      };    
      
      /** Constructor */
      RungeKuttaEngine(const Config& rkConfig);

      virtual ~RungeKuttaEngine();

      /** resolve the boundary situation - for charged particles */
      ExtrapolationCode propagate(ExCellCharged& ecCell,
                                  const Surface& sf,
                                  PropDirection dir=alongMomentum,
                                  const BoundaryCheck& bcheck = true,
                                  bool returnCurvilinear = true) const final;  

      /** resolve the boundary situation - for neutral particles */
      ExtrapolationCode propagate(ExCellNeutral& enCell,
                                  const Surface& sf,
                                  PropDirection dir=alongMomentum,
                                  const BoundaryCheck& bcheck = true,
                                  bool returnCurvilinear = true) const final;
                                  
      /** Set configuration method */
      void setConfiguration(const Config& meConfig);
      
      /** Get configuration method */
      Config getConfiguration() const;                                    
    
    protected:
      Config            m_config;  //!< configuration class
      
      RungeKuttaUtils   m_rkUtils; //!< RungeKuttaUtils class
      
    private:
      const Logger& logger() const {return *m_config.logger;}

    /** Templated RungeKutta propagation method - charged/neutral */
      template <class T> bool propagateRungeKuttaT(ExtrapolationCell<T>& eCell,
                                                   PropagationCache& pCache,
                                                   const T& tParameters,
                                                   const Surface& sf) const; 
        
      /** Internal RungeKutta propagation method for propation with jacobian*/     
      bool propagateWithJacobian(int navigationStep,
                                 PropagationCache& pCache,
                                 int surfaceType,
	                             double* sVector) const;

      /** Propagation methods runge kutta step - returns the step length */
      double rungeKuttaStep(int navigationStep,
                            PropagationCache& pCache,
                            double,
                            bool&) const;

      /** Propagation methods runge kutta step - returns the step length*/
      double rungeKuttaStepWithGradient(int navigationStep,    
                                        PropagationCache& pCache,
                                        double,
                                        bool& ) const;

      /** Propagation methods straight line step*/
      double straightLineStep(int navigationStep, PropagationCache& pCache, double) const;

      /** Step estimator with directions correction */      
      double stepEstimatorWithCurvature(PropagationCache& pCache, int, double*, bool&) const;

      /** Build new track parameters without propagation */
      const TrackParameters* buildTrackParametersWithoutPropagation(const TrackParameters &, double*) const;

      /** Build new track parameters without propagation */
      const NeutralParameters* buildNeutralParametersWithoutPropagation(const NeutralParameters&, double*) const;

      /** Test new propagation to cylinder boundary */
      bool newCrossPoint(const CylinderSurface&, const double *, const double  *) const;
       
      // get the field - with the fast option 
      void getField        (double*,double*        ) const;
      void getFieldGradient(double*,double*,double*) const;
      
   };

  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods for magnetic field information
  /////////////////////////////////////////////////////////////////////////////////

  inline void RungeKuttaEngine::getField(double* R,double* H) const
  {
    m_config.fieldService->getField(R,H);
  }

  inline void RungeKuttaEngine::getFieldGradient(double* R,double* H,double* dH) const
  {
    m_config.fieldService->getField(R,H,dH);
  }

  /** Return the configuration object */    
  inline RungeKuttaEngine::Config RungeKuttaEngine::getConfiguration() const { return m_config; }

  ////////////////////////////////////////////////////////////////////////////////
  // Templated method
  ////////////////////////////////////////////////////////////////////////////////
#include "ACTS/Extrapolation/detail/RungeKuttaEngine.icc"
  
}
 
#endif // ACTS_RUNGEKUTTAENGINE_RUNGEKUTAENGINE_H
