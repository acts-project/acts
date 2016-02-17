/////////////////////////////////////////////////////////////////////////////////
//  Header file for class RungeKuttaEngine, ATS project
/////////////////////////////////////////////////////////////////////////////////

#ifndef ATS_RUNGEKUTTAENGINE_RUNGEKUTAENGINE_H
#define ATS_RUNGEKUTTAENGINE_RUNGEKUTAENGINE_H 1

// Gaudi
#include "GaudiKernel/ServiceHandle.h"
// Core module
#include "CoreInterfaces/ServiceBase.h"
// Extrapolation module
#include "ExtrapolationInterfaces/IPropagationEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "ExtrapolationUtils/RungeKuttaUtils.h"
// EventData module
#include "EventDataUtils/PropDirection.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"
// Geometry module
#include "Surfaces/Surface.h"
#include "Surfaces/CylinderSurface.h"
#include "Surfaces/ConeSurface.h"
#include "Surfaces/BoundaryCheck.h"
// MagneticField module
#include "MagneticFieldUtils/MagneticFieldProperties.h"
#include "MagneticFieldInterfaces/IMagneticFieldSvc.h"

namespace Ats {

  /**
    @struct PropagationCache 
    
     Helper struct to allow state-less propagation.
    
    @author Andreas.Salzburger -at- cern.ch 
    
    */
    
  struct PropagationCache { 
      
      // configuration
      double                    direction;
      BoundaryCheck             boundaryCheck;
      MagneticFieldProperties   mFieldMode;
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
      AtsSymMatrixD<5>*          covariance; 
      double                    jacobian[25];
      
      
      PropagationCache() :
       direction(alongMomentum),
       boundaryCheck(true),
       mFieldMode(FullField),
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

  class RungeKuttaEngine : public ServiceBase, virtual public IPropagationEngine {
      
    public:
      
      /** Constructor */
      RungeKuttaEngine(const std::string& name, ISvcLocator* svc);

      virtual ~RungeKuttaEngine();

      /** Destructor */
      StatusCode initialize() final;

      /** AlgTool finalize method */
      StatusCode finalize() final;
      
      /** Query the interfaces **/
      StatusCode queryInterface( const InterfaceID& riid, void** ppvInterface );

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
                                  
      
    private:
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
      
      /////////////////////////////////////////////////////////////////////////////////
      // Private data members: 
      /////////////////////////////////////////////////////////////////////////////////
      ServiceHandle<IMagneticFieldSvc>  m_fieldService;  //!< the field service 

      double m_dlt                                    ;  //!< accuracy parameter
      double m_helixStep                              ;  //!< max step whith helix model
      double m_straightStep                           ;  //!< max step whith srtaight line model
      double m_maxPathLength                          ;  //!< max overal path length 
      bool   m_usegradient                            ;  //!< use magnetif field gradient
      
   };

  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods for magnetic field information
  /////////////////////////////////////////////////////////////////////////////////

  inline void RungeKuttaEngine::getField(double* R,double* H) const
  {
    m_fieldService->getField(R,H);
  }

  inline void RungeKuttaEngine::getFieldGradient(double* R,double* H,double* dH) const
  {
    m_fieldService->getField(R,H,dH);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Templated method
  ////////////////////////////////////////////////////////////////////////////////
#include "RungeKuttaEngine/RungeKuttaEngine.icc"
  
}
 
#endif // ATS_RUNGEKUTTAENGINE_RUNGEKUTAENGINE_H
