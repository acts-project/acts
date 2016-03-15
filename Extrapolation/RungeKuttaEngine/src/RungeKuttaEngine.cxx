/////////////////////////////////////////////////////////////////////////////////
// RungeKuttaEngine.cxx, ACTS project
/////////////////////////////////////////////////////////////////////////////////

// Event module
#include "EventDataUtils/CoordinateTransformations.h"
// Extrapolation module
#include "RungeKuttaEngine/RungeKuttaEngine.h"
#include "ExtrapolationInterfaces/ExtrapolationMacros.h"
#include "ExtrapolationUtils/TransportJacobian.h"
// Geometry module
#include "Surfaces/Surface.h"
#include "Surfaces/DiscSurface.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/PerigeeSurface.h"
#include "Surfaces/StraightLineSurface.h"
// MagneticField module
#include "MagneticFieldUtils/MagneticFieldProperties.h"

/////////////////////////////////////////////////////////////////////////////////
// Constructor
/////////////////////////////////////////////////////////////////////////////////
Acts::RungeKuttaEngine::RungeKuttaEngine(const std::string& name, ISvcLocator* svc) :  
  ServiceBase(name, svc),
  m_fieldService("", name),
  m_dlt(0.000200),
  m_helixStep(1.), 
  m_straightStep(.01),
  m_maxPathLength(25000.),
  m_usegradient(false)
  
{
  // steering of the screen outoput (SOP)
  declareProperty("OutputPrefix"       , m_sopPrefix);
  declareProperty("OutputPostfix"      , m_sopPostfix);
  // property declaration
  declareProperty("AccuracyParameter"  , m_dlt          );
  declareProperty("MaxHelixStep"       , m_helixStep    );
  declareProperty("MaxStraightLineStep", m_straightStep);
  declareProperty("MaxPathLength"      , m_maxPathLength);
  declareProperty("IncludeBgradients"  , m_usegradient );
  declareProperty("MagneticFieldSvc"   , m_fieldService);
}

/////////////////////////////////////////////////////////////////////////////////
// Destructor
/////////////////////////////////////////////////////////////////////////////////
Acts::RungeKuttaEngine::~RungeKuttaEngine(){}

/////////////////////////////////////////////////////////////////////////////////
// Query the interface
/////////////////////////////////////////////////////////////////////////////////
StatusCode Acts::RungeKuttaEngine::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
  if ( IID_IPropagationEngine == riid )
    *ppvInterface = (IPropagationEngine*)this;
  else  {
    // Interface is not directly available: try out a base class
    return Service::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}  

/////////////////////////////////////////////////////////////////////////////////
// initialize
/////////////////////////////////////////////////////////////////////////////////
StatusCode Acts::RungeKuttaEngine::initialize()
{
  MSG_DEBUG("initialize()");
  // retrieve the tracking geometry servcie - crucial, abort when it can not be retrieved
  RETRIEVE_FATAL(m_fieldService);
  return StatusCode::SUCCESS;
}

/////////////////////////////////////////////////////////////////////////////////
// finalize
/////////////////////////////////////////////////////////////////////////////////
StatusCode Acts::RungeKuttaEngine::finalize()
{
  MSG_DEBUG("finalize()");
  return StatusCode::SUCCESS;
}

/////////////////////////////////////////////////////////////////////////////////
// Main function for NeutralParameters propagation 
/////////////////////////////////////////////////////////////////////////////////
Acts::ExtrapolationCode Acts::RungeKuttaEngine::propagate(ExCellNeutral& eCell,
                                                        const Surface& sf,
                                                        PropDirection pDir,
                                                        const BoundaryCheck& bcheck,
                                                        bool returnCurvilinear) const
{ 
    EX_MSG_DEBUG(++eCell.navigationStep, "propagate", "neut", "propagation engine called with neutral parameters with propagation direction " << pDir );     
    // it is the final propagation if it is the endSurface
    bool finalPropagation = (eCell.endSurface == (&sf));

    // create the PropagationCache 
    PropagationCache pCache;
    
    // get the start parameters
    const NeutralParameters* sParameters = eCell.leadParameters;
    const NeutralParameters* nParameters = nullptr;

    // if the desination surface is the start surface -> bail out and build parameters directly
    if (&sf == &(sParameters->associatedSurface())){ 
       EX_MSG_VERBOSE(eCell.navigationStep, "propagate", "neut", "parameters are already on the surface, returning.");             
       nParameters = buildNeutralParametersWithoutPropagation(*sParameters,pCache.jacobian);
       // fast exit
       eCell.lastLeadParameters = sParameters;
       // assign the lead and end parameters
       eCell.leadParameters = nParameters;
       // check if the propagation was called with directly, then lead parameters become end parameters
       if (eCell.checkConfigurationMode(ExtrapolationMode::Direct)) 
           eCell.endParameters = eCell.leadParameters;
       // return success or in progress
       return (finalPropagation ? ExtrapolationCode::SuccessDestination : ExtrapolationCode::InProgress);
    }  
    // specify the parameters for the propagation    
    pCache.maxPathLength     = eCell.pathLimit < 0. ? m_maxPathLength : (eCell.pathLimit - eCell.pathLength);
    pCache.direction         = double(pDir);
    pCache.boundaryCheck     = bcheck;
    pCache.returnCurvilinear = returnCurvilinear;
    pCache.useJacobian       = eCell.leadParameters->covariance();
    // neutral transport sets mconditions to false
    pCache.mcondition        = false;

    // the result through propagation
    if (propagateRungeKuttaT<NeutralParameters>(eCell, pCache, *sParameters, sf)){
        // screen output
        EX_MSG_VERBOSE(eCell.navigationStep, "propagate", "neut", "propagation to surface was successful.");     
        // create a new covariance matrix
        std::unique_ptr<ActsSymMatrixD<5> > cov;
        if(pCache.covariance)
          cov.reset(new ActsSymMatrixD<5>(*pCache.covariance));
        
        // create the new parameters
        if (!pCache.returnCurvilinear){
          // new parameters bound to the surface
          ActsVectorD<5> pars;
          pars << pCache.parameters[0],pCache.parameters[1],pCache.parameters[2],pCache.parameters[3],pCache.parameters[4];
          nParameters = new NeutralBoundParameters(std::move(cov),std::move(pars),sf);
        } else {
          // new curvilinear parameters
          Acts::Vector3D gp(pCache.pVector[0],pCache.pVector[1],pCache.pVector[2]);
          Acts::Vector3D mom(pCache.pVector[3],pCache.pVector[4],pCache.pVector[5]);
          mom /= fabs(pCache.parameters[4]);
          nParameters = new NeutralCurvilinearParameters(std::move(cov),gp,mom);
        }
    } 
    // only go on if the parameter creation worked 
    if (nParameters){        
        // @!TODO fill the jacobian
        // if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::CollectJacobians)){ doSomething; }
        
        
        // fill the transport information - only if the propation direction is not 0 ('anyDirection')
        if (pDir!=anyDirection){
           EX_MSG_VERBOSE(eCell.navigationStep,"propagate", "neut", "path length of " << pCache.step << " added to the extrapolation cell (limit = " << eCell.pathLimit << ")" );    
           eCell.stepTransport(sf,pCache.step);
        }
		// cache the last lead parameters
		eCell.lastLeadParameters = sParameters;
        // now exchange the lead parameters 
        // create the new curvilinear paramters at the surface intersection -> if so, trigger the success
        eCell.leadParameters = nParameters;
        
        // now check if it is valid it's further away than the pathLimit
        if (eCell.pathLimitReached(m_dlt)){
            // screen output
            EX_MSG_VERBOSE(eCell.navigationStep,"propagate", "neut", "path limit of " << eCell.pathLimit << " reached. Stopping extrapolation."); 
            return ExtrapolationCode::SuccessPathLimit;
        }  
                                                       
        // check if the propagation was called with directly, then lead parameters become end parameters
        if (eCell.checkConfigurationMode(ExtrapolationMode::Direct)) 
	        eCell.endParameters = eCell.leadParameters;
	    // return success for the final destination or in progress                                                                   
        return (finalPropagation ? ExtrapolationCode::SuccessDestination : ExtrapolationCode::InProgress);
    } else {
        // give some screen output
        EX_MSG_VERBOSE(eCell.navigationStep,"propagate", "neut", "intersection with the surface did not succeed.");
    }                                                                     
   // return - recovered means that the leadParameters are the input ones 
   return (finalPropagation ? ExtrapolationCode::FailureDestination : ExtrapolationCode::Recovered) ;
}

/////////////////////////////////////////////////////////////////////////////////
// Main function for TrackParameters propagation 
/////////////////////////////////////////////////////////////////////////////////
Acts::ExtrapolationCode Acts::RungeKuttaEngine::propagate(ExCellCharged& eCell,
                                                        const Surface& sf,
                                                        PropDirection pDir,
                                                        const BoundaryCheck& bcheck,
                                                        bool returnCurvilinear) const
{
    EX_MSG_DEBUG(++eCell.navigationStep, "propagate", "char", "propagation engine called with charged parameters with propagation direction " << pDir ); 
    // it is the final propagation if it is the endSurface
    bool finalPropagation = (eCell.endSurface == (&sf));
    
    // the start and teh result
    const TrackParameters* pParameters = nullptr;
    const TrackParameters* sParameters = eCell.leadParameters;
    
    // build the propagation cache
    PropagationCache pCache;
    
    // if the desination surface is the start surface -> bail out and build parameters directly
    if (&sf == &(eCell.leadParameters->associatedSurface())) { 
        EX_MSG_VERBOSE(eCell.navigationStep, "propagate", "char", "parameters are already on the surface, returning.");             
       // create them without propagation -> success     
       pParameters = buildTrackParametersWithoutPropagation(*sParameters,pCache.jacobian);
       // fast exit
       eCell.lastLeadParameters = sParameters;
       // assign the lead and end parameters
       eCell.leadParameters = pParameters;
       // check if the propagation was called with directly, then lead parameters become end parameters
       if (eCell.checkConfigurationMode(ExtrapolationMode::Direct)) 
           eCell.endParameters = eCell.leadParameters;
       // return success or in progress
       return (finalPropagation ? ExtrapolationCode::SuccessDestination : ExtrapolationCode::InProgress);
    }

    // and configure the propagation cache now
    pCache.maxPathLength     = eCell.pathLimit < 0. ? m_maxPathLength : (eCell.pathLimit - eCell.pathLength);
    pCache.direction         = double(pDir);
    pCache.boundaryCheck     = bcheck;
    pCache.returnCurvilinear = returnCurvilinear;
    pCache.useJacobian       = eCell.leadParameters->covariance(); 
    pCache.mFieldMode        = eCell.mFieldMode;
    pCache.mcondition        = (eCell.mFieldMode.magneticFieldMode() != 0 ) ? true : false;
    pCache.needgradient      = (pCache.useJacobian && m_usegradient ) ? true : false;
    
    // propagate with templated helper function
    if (propagateRungeKuttaT<TrackParameters>(eCell, pCache, *sParameters, sf)){
        EX_MSG_VERBOSE(eCell.navigationStep, "propagate", "char", "propagation to surface was successful.");     
        // create the new parameters
        // create a new covariance matrix
        std::unique_ptr<ActsSymMatrixD<5> > cov;
        if(pCache.covariance)
          cov.reset(new ActsSymMatrixD<5>(*pCache.covariance));
        // create the parameter vector
        ActsVectorD<5> pars;
        pars << pCache.parameters[0],pCache.parameters[1],pCache.parameters[2],pCache.parameters[3],pCache.parameters[4];
        // create the new parameters
        if (!pCache.returnCurvilinear){
          // new parameters bound to the surface
          pParameters = new BoundParameters(std::move(cov),std::move(pars),sf);
        } else {
          // get the charge
          double charge = pCache.parameters[4] > 0. ? 1. : -1.;
          // new curvilinear parameters
          Acts::Vector3D gp(pCache.pVector[0],pCache.pVector[1],pCache.pVector[2]);       
          pParameters = new CurvilinearParameters(std::move(cov),gp,coordinate_transformation::parameters2globalMomentum(pars),charge);
        }
    }
    // set the return type according to how the propagation went
    if (pParameters){
       // @!TODO fill the jacobian
        // if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::CollectJacobians)){ doSomething; }
        
       // cache the last lead parameters, useful in case a navigation error occured
       eCell.lastLeadParameters = eCell.leadParameters;
       // assign the lead and end parameters
       eCell.leadParameters = pParameters;
       // check what to do with the path Length
       if (eCell.checkConfigurationMode(ExtrapolationMode::StopWithPathLimit) || eCell.pathLength > 0){
           // add the new propagation length to the path length
           eCell.pathLength += pCache.step;
           // check if Limit reached
           if (eCell.pathLimitReached(m_dlt)){
               EX_MSG_VERBOSE(eCell.navigationStep, "propagate", "char", "path limit of " << eCell.pathLimit << " successfully reached -> stopping." ); 
               return ExtrapolationCode::SuccessPathLimit;
           }
       }
       // check if the propagation was called with directly, then lead parameters become end parameters
       if (eCell.checkConfigurationMode(ExtrapolationMode::Direct)) 
           eCell.endParameters = eCell.leadParameters;
	 
       // return Success only if it is the final propagation - the extrapolation engine knows that 
       return (finalPropagation ? ExtrapolationCode::SuccessDestination : ExtrapolationCode::InProgress);
   }                                                                      
   // return - recovered means that the leadParameters are the input ones 
   return (finalPropagation ? Acts::ExtrapolationCode::FailureDestination : Acts::ExtrapolationCode::Recovered) ;
}
                                              
/////////////////////////////////////////////////////////////////////////////////
// Runge Kutta main program for propagation with or without Jacobian
/////////////////////////////////////////////////////////////////////////////////
bool Acts::RungeKuttaEngine::propagateWithJacobian(int navigationStep, PropagationCache& pCache, int kind, double* Su) const
{
  
  EX_MSG_VERBOSE(navigationStep, "propagate", "<T> ", "propagateWithJacobian called with  internal surface type " << kind ); 
  
  const double Smax        = 1000.                ;  // max. step allowed
  double       Wwrong      = 500.                 ;  // Max way with wrong direction
  double*      R           = &(pCache.pVector[ 0]);  // Start coordinates
  double*      A           = &(pCache.pVector[ 3]);  // Start directions
  double*      SA          = &(pCache.pVector[42]); SA[0]=SA[1]=SA[2]=0.;
  pCache.maxPathLimit      = false           ;  

  if (pCache.mcondition && fabs(pCache.pVector[6]) > .1) return false; 

  // Step estimation until surface
  //
  Acts::RungeKuttaUtils utils;
  bool Q; double S, step = utils.stepEstimator(kind,Su,pCache.pVector,Q); if(!Q) return false; 

  bool dir = true;
  if (pCache.mcondition && pCache.direction && pCache.direction*step < 0.)  {
    step = -step; dir = false;
  }

  step>Smax ? S=Smax : step<-Smax ? S=-Smax : S=step;
  double So = fabs(S); int iS = 0;

  bool InS = false;

  // Rkuta extrapolation
  //
  int niter = 0;
  pCache.newfield = true;
  
  // whie loop over the steps
  while (fabs(step) > m_straightStep) {

    // maximum number of steps
    if (++niter > 10000) {
        //!< @TODO make max number configurable
        EX_MSG_DEBUG(navigationStep, "propagate", "<T> ", "maximimum number of integration steps (" << 10000 <<") reached. Aborting."  ); 
        return false; 
    }

    // propagation in magnetic field  - with or without field ( or gradient )
    if(pCache.mcondition) {
       if(!pCache.needgradient) pCache.step += (S=rungeKuttaStep            (navigationStep,pCache,S,InS));
       else                     pCache.step += (S=rungeKuttaStepWithGradient(navigationStep,pCache,S,InS));
    }
    else  { // or within straight line 
      pCache.step += (S=straightLineStep(navigationStep,pCache,S));
    }

    step = stepEstimatorWithCurvature(pCache,kind,Su,Q); 
    
    if (!Q) {
        EX_MSG_DEBUG(navigationStep, "propagate", "<T> ", "step estimation with curvature did not succeed. Aborting"  ); 
        return false; 
    }

    if(!dir) {
      if (pCache.direction && pCache.direction*pCache.step < 0.)  step = -step;
      else dir  =  true;
    }
    
    if (S*step<0.) {S = -S; ++iS;}
 
    // check if the step made sense 
    double aS    = fabs(S   );
    double aStep = fabs(step);
    if     (    aS > aStep             )  S = step;
    else if(!iS && InS && aS*2. < aStep)  S*=2.   ;
    
    if (!dir && fabs(pCache.step) > Wwrong ) {
        EX_MSG_DEBUG(navigationStep, "propagate", "<T> ", "step into the wrong direction done. Aborting"  ); 
        return false; 
    }
    
    if(iS > 10 || (iS>3 && fabs(S)>=So)) { 
        if(!kind) break; 
        EX_MSG_DEBUG(navigationStep, "propagate", "<T> ", "Abort triggered."  ); 
        return false;
    
    }
        
    double dW =  pCache.maxPathLength-fabs(pCache.step);
    if(fabs(S) > dW) {S > 0. ? S = dW : S = -dW; step = S; pCache.maxPathLimit = true;}

    So=fabs(S);
        
  } // end of while loop
  
  EX_MSG_VERBOSE(navigationStep, "propagate", "<T> ", "numerical integration is done." ); 
  
  // Output track parameteres
  //
  pCache.step += step;

  if (fabs(step) < .001) return true;

  A [0]+=(SA[0]*step); 
  A [1]+=(SA[1]*step);
  A [2]+=(SA[2]*step);
  double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);

  R[0]+=step*(A[0]-.5*step*SA[0]); A[0]*=CBA;
  R[1]+=step*(A[1]-.5*step*SA[1]); A[1]*=CBA;
  R[2]+=step*(A[2]-.5*step*SA[2]); A[2]*=CBA;
  
 
  return true;
}

/////////////////////////////////////////////////////////////////////////////////
// Runge Kutta trajectory model (units->mm,MeV,kGauss)
// Uses Nystroem algorithm (See Handbook Net. Bur. ofStandards, procedure 25.5.20)
/////////////////////////////////////////////////////////////////////////////////
double Acts::RungeKuttaEngine::rungeKuttaStep(int navigationStep, PropagationCache& pCache, double S,  bool& InS) const
{
  
  EX_MSG_VERBOSE(navigationStep, "propagate", "<T> ", "rungeKuttaStep called"); 
  
  bool Jac = pCache.useJacobian;
    
  double* R    =          &(pCache.pVector[ 0]);            // Coordinates 
  double* A    =          &(pCache.pVector[ 3]);            // Directions
  double* sA   =          &(pCache.pVector[42]);
  double  Pi   =  149.89626*pCache.pVector[6];            // Invert mometum/2. 
  double  dltm = m_dlt*.03      ;

  double f0[3],f[3]; 
  
  // if new field is required get it
  if (pCache.newfield) getField(R,f0); 
  else { f0[0]=pCache.field[0]; f0[1]=pCache.field[1]; f0[2]=pCache.field[2];}

  bool Helix = false; if (fabs(S) < m_helixStep) Helix = true; 
  
  while(S != 0.) {
     
    double S3=(1./3.)*S, S4=.25*S, PS2=Pi*S;

    // First point
    //   
    double H0[3] = {f0[0]*PS2, f0[1]*PS2, f0[2]*PS2};
    double A0    = A[1]*H0[2]-A[2]*H0[1]            ;
    double B0    = A[2]*H0[0]-A[0]*H0[2]            ;
    double C0    = A[0]*H0[1]-A[1]*H0[0]            ;
    double A2    = A0+A[0]                          ;
    double B2    = B0+A[1]                          ;
    double C2    = C0+A[2]                          ;
    double A1    = A2+A[0]                          ;
    double B1    = B2+A[1]                          ;
    double C1    = C2+A[2]                          ;
    
    // Second point
    //
    if(!Helix) {
      double gP[3]={R[0]+A1*S4, R[1]+B1*S4, R[2]+C1*S4};
      getField(gP,f);
    }
    else {f[0]=f0[0]; f[1]=f0[1]; f[2]=f0[2];}

    double H1[3] = {f[0]*PS2,f[1]*PS2,f[2]*PS2}; 
    double A3    = (A[0]+B2*H1[2])-C2*H1[1]    ; 
    double B3    = (A[1]+C2*H1[0])-A2*H1[2]    ; 
    double C3    = (A[2]+A2*H1[1])-B2*H1[0]    ;
    double A4    = (A[0]+B3*H1[2])-C3*H1[1]    ; 
    double B4    = (A[1]+C3*H1[0])-A3*H1[2]    ; 
    double C4    = (A[2]+A3*H1[1])-B3*H1[0]    ;
    double A5    = 2.*A4-A[0]                  ; 
    double B5    = 2.*B4-A[1]                  ; 
    double C5    = 2.*C4-A[2]                  ;    

    // Last point
    //
    if(!Helix) {
      double gP[3]={R[0]+S*A4, R[1]+S*B4, R[2]+S*C4};    
      getField(gP,f);
    }
    else       {f[0]=f0[0]; f[1]=f0[1]; f[2]=f0[2];} 

    double H2[3] = {f[0]*PS2,f[1]*PS2,f[2]*PS2}; 
    double A6    = B5*H2[2]-C5*H2[1]           ;
    double B6    = C5*H2[0]-A5*H2[2]           ;
    double C6    = A5*H2[1]-B5*H2[0]           ;
    
    // Test approximation quality on give step and possible step reduction
    //
    double EST = fabs((A1+A6)-(A3+A4))+fabs((B1+B6)-(B3+B4))+fabs((C1+C6)-(C3+C4)); 
    if(EST>m_dlt) {S*=.5; dltm = 0.; continue;} EST<dltm ? InS = true : InS = false;

    // Parameters calculation
    //   
    double A00 = A[0], A11=A[1], A22=A[2];

    A[0] = 2.*A3+(A0+A5+A6); 
    A[1] = 2.*B3+(B0+B5+B6); 
    A[2] = 2.*C3+(C0+C5+C6);
    
    double D  = (A[0]*A[0]+A[1]*A[1])+(A[2]*A[2]-9.);
    double Sl = 2./S                                ;
    D         = (1./3.)-((1./648.)*D)*(12.-D)       ;

    R[0] +=(A2+A3+A4)*S3;
    R[1] +=(B2+B3+B4)*S3;
    R[2] +=(C2+C3+C4)*S3;
    A[0] *=D            ;
    A[1] *=D            ;
    A[2] *=D            ;
    sA[0] = A6*Sl       ; 
    sA[1] = B6*Sl       ;
    sA[2] = C6*Sl       ; 

    pCache.field[0]=f[0]; pCache.field[1]=f[1]; pCache.field[2]=f[2]; pCache.newfield = false;

    if(!Jac) return S;

    // Jacobian calculation
    //
    double* d2A = &pCache.pVector[24];
    double* d3A = &pCache.pVector[31]; 
    double* d4A = &pCache.pVector[38]; 
    double d2A0 = H0[2]*d2A[1]-H0[1]*d2A[2];
    double d2B0 = H0[0]*d2A[2]-H0[2]*d2A[0];
    double d2C0 = H0[1]*d2A[0]-H0[0]*d2A[1];
    double d3A0 = H0[2]*d3A[1]-H0[1]*d3A[2];
    double d3B0 = H0[0]*d3A[2]-H0[2]*d3A[0];
    double d3C0 = H0[1]*d3A[0]-H0[0]*d3A[1];
    double d4A0 =(A0+H0[2]*d4A[1])-H0[1]*d4A[2];
    double d4B0 =(B0+H0[0]*d4A[2])-H0[2]*d4A[0];
    double d4C0 =(C0+H0[1]*d4A[0])-H0[0]*d4A[1];
    double d2A2 = d2A0+d2A[0];                
    double d2B2 = d2B0+d2A[1];                
    double d2C2 = d2C0+d2A[2];
    double d3A2 = d3A0+d3A[0];                
    double d3B2 = d3B0+d3A[1];                
    double d3C2 = d3C0+d3A[2];
    double d4A2 = d4A0+d4A[0];                
    double d4B2 = d4B0+d4A[1];                
    double d4C2 = d4C0+d4A[2];
    double d0   = d4A[0]-A00;
    double d1   = d4A[1]-A11;
    double d2   = d4A[2]-A22;
    double d2A3 = ( d2A[0]+d2B2*H1[2])-d2C2*H1[1];
    double d2B3 = ( d2A[1]+d2C2*H1[0])-d2A2*H1[2];
    double d2C3 = ( d2A[2]+d2A2*H1[1])-d2B2*H1[0];
    double d3A3 = ( d3A[0]+d3B2*H1[2])-d3C2*H1[1];
    double d3B3 = ( d3A[1]+d3C2*H1[0])-d3A2*H1[2];
    double d3C3 = ( d3A[2]+d3A2*H1[1])-d3B2*H1[0];
    double d4A3 = ((A3+d0)+d4B2*H1[2])-d4C2*H1[1];
    double d4B3 = ((B3+d1)+d4C2*H1[0])-d4A2*H1[2];
    double d4C3 = ((C3+d2)+d4A2*H1[1])-d4B2*H1[0];
    double d2A4 = ( d2A[0]+d2B3*H1[2])-d2C3*H1[1];
    double d2B4 = ( d2A[1]+d2C3*H1[0])-d2A3*H1[2];
    double d2C4 = ( d2A[2]+d2A3*H1[1])-d2B3*H1[0];
    double d3A4 = ( d3A[0]+d3B3*H1[2])-d3C3*H1[1];
    double d3B4 = ( d3A[1]+d3C3*H1[0])-d3A3*H1[2];
    double d3C4 = ( d3A[2]+d3A3*H1[1])-d3B3*H1[0];
    double d4A4 = ((A4+d0)+d4B3*H1[2])-d4C3*H1[1];
    double d4B4 = ((B4+d1)+d4C3*H1[0])-d4A3*H1[2];
    double d4C4 = ((C4+d2)+d4A3*H1[1])-d4B3*H1[0];
    double d2A5 = 2.*d2A4-d2A[0];            
    double d2B5 = 2.*d2B4-d2A[1];            
    double d2C5 = 2.*d2C4-d2A[2];
    double d3A5 = 2.*d3A4-d3A[0];            
    double d3B5 = 2.*d3B4-d3A[1];            
    double d3C5 = 2.*d3C4-d3A[2];            
    double d4A5 = 2.*d4A4-d4A[0];            
    double d4B5 = 2.*d4B4-d4A[1];            
    double d4C5 = 2.*d4C4-d4A[2];            
    double d2A6 = d2B5*H2[2]-d2C5*H2[1];      
    double d2B6 = d2C5*H2[0]-d2A5*H2[2];      
    double d2C6 = d2A5*H2[1]-d2B5*H2[0];      
    double d3A6 = d3B5*H2[2]-d3C5*H2[1];      
    double d3B6 = d3C5*H2[0]-d3A5*H2[2];      
    double d3C6 = d3A5*H2[1]-d3B5*H2[0];
    double d4A6 = d4B5*H2[2]-d4C5*H2[1];      
    double d4B6 = d4C5*H2[0]-d4A5*H2[2];      
    double d4C6 = d4A5*H2[1]-d4B5*H2[0];      
      
    double* dR  = &pCache.pVector[21];
    dR [0]+=(d2A2+d2A3+d2A4)*S3;
    dR [1]+=(d2B2+d2B3+d2B4)*S3;
    dR [2]+=(d2C2+d2C3+d2C4)*S3;
    d2A[0] =((d2A0+2.*d2A3)+(d2A5+d2A6))*(1./3.);      
    d2A[1] =((d2B0+2.*d2B3)+(d2B5+d2B6))*(1./3.); 
    d2A[2] =((d2C0+2.*d2C3)+(d2C5+d2C6))*(1./3.);

    dR          = &pCache.pVector[28];
    dR [0]+=(d3A2+d3A3+d3A4)*S3;
    dR [1]+=(d3B2+d3B3+d3B4)*S3;
    dR [2]+=(d3C2+d3C3+d3C4)*S3;
    d3A[0] =((d3A0+2.*d3A3)+(d3A5+d3A6))*(1./3.);      
    d3A[1] =((d3B0+2.*d3B3)+(d3B5+d3B6))*(1./3.); 
    d3A[2] =((d3C0+2.*d3C3)+(d3C5+d3C6))*(1./3.);

    dR          = &pCache.pVector[35];
    dR [0]+=(d4A2+d4A3+d4A4)*S3;
    dR [1]+=(d4B2+d4B3+d4B4)*S3;
    dR [2]+=(d4C2+d4C3+d4C4)*S3;
    d4A[0] =((d4A0+2.*d4A3)+(d4A5+d4A6+A6))*(1./3.);      
    d4A[1] =((d4B0+2.*d4B3)+(d4B5+d4B6+B6))*(1./3.); 
    d4A[2] =((d4C0+2.*d4C3)+(d4C5+d4C6+C6))*(1./3.);
    return S;
  }
  return S;
}

/////////////////////////////////////////////////////////////////////////////////
// Runge Kutta trajectory model (units->mm,MeV,kGauss)
// Uses Nystroem algorithm (See Handbook Net. Bur. ofStandards, procedure 25.5.20)
//    Where magnetic field information iS              
//    f[ 0],f[ 1],f[ 2] - Hx    , Hy     and Hz of the magnetic field         
//    f[ 3],f[ 4],f[ 5] - dHx/dx, dHx/dy and dHx/dz                           
//    f[ 6],f[ 7],f[ 8] - dHy/dx, dHy/dy and dHy/dz                           
//    f[ 9],f[10],f[11] - dHz/dx, dHz/dy and dHz/dz                           
//                                                                                   
/////////////////////////////////////////////////////////////////////////////////
double Acts::RungeKuttaEngine::rungeKuttaStepWithGradient(int navigationStep, PropagationCache& pCache, double S, bool& InS) const
{
  
  EX_MSG_VERBOSE(navigationStep, "propagate", "<T> ", "rungeKuttaStepWithGradient called"); 
  
  const double C33 = 1./3.      ;  
  double* R    =          &(pCache.pVector[ 0]);           // Coordinates 
  double* A    =          &(pCache.pVector[ 3]);           // Directions
  double* sA   =          &(pCache.pVector[42]);
  double  Pi   =  149.89626*pCache.pVector[6];           // Invert mometum/2. 
  double  dltm = m_dlt*.03      ;

  double f0[3],f1[3],f2[3],g0[9],g1[9],g2[9],H0[12],H1[12],H2[12];
  getFieldGradient(R,f0,g0);

  while (S != 0.) {
 
    
    double S3=C33*S, S4=.25*S, PS2=Pi*S;

    // First point
    //   
    H0[0] = f0[0]*PS2; H0[1] = f0[1]*PS2; H0[2] = f0[2]*PS2;
    double A0    = A[1]*H0[2]-A[2]*H0[1]             ;
    double B0    = A[2]*H0[0]-A[0]*H0[2]             ;
    double C0    = A[0]*H0[1]-A[1]*H0[0]             ;
    double A2    = A[0]+A0                           ;
    double B2    = A[1]+B0                           ;
    double C2    = A[2]+C0                           ;
    double A1    = A2+A[0]                           ;
    double B1    = B2+A[1]                           ;
    double C1    = C2+A[2]                           ;
    
    // Second point
    //
    double gP1[3]={R[0]+A1*S4, R[1]+B1*S4, R[2]+C1*S4};
    getFieldGradient(gP1,f1,g1);

    H1[0] = f1[0]*PS2; H1[1] = f1[1]*PS2; H1[2] = f1[2]*PS2; 
    double A3    = B2*H1[2]-C2*H1[1]+A[0]         ; 
    double B3    = C2*H1[0]-A2*H1[2]+A[1]         ; 
    double C3    = A2*H1[1]-B2*H1[0]+A[2]         ;
    double A4    = B3*H1[2]-C3*H1[1]+A[0]         ; 
    double B4    = C3*H1[0]-A3*H1[2]+A[1]         ; 
    double C4    = A3*H1[1]-B3*H1[0]+A[2]         ;
    double A5    = A4-A[0]+A4                     ; 
    double B5    = B4-A[1]+B4                     ; 
    double C5    = C4-A[2]+C4                     ;
    
    // Last point
    //
    double gP2[3]={R[0]+S*A4, R[1]+S*B4, R[2]+S*C4};    
    getFieldGradient(gP2,f2,g2);

    H2[0] = f2[0]*PS2; H2[1] = f2[1]*PS2; H2[2] = f2[2]*PS2; 
    double A6    = B5*H2[2]-C5*H2[1]              ;
    double B6    = C5*H2[0]-A5*H2[2]              ;
    double C6    = A5*H2[1]-B5*H2[0]              ;
    
    // Test approximation quality on give step and possible step reduction
    //
    double EST = fabs((A1+A6)-(A3+A4))+fabs((B1+B6)-(B3+B4))+fabs((C1+C6)-(C3+C4)); 
    if(EST>m_dlt) {S*=.5; dltm = 0.; continue;} EST<dltm ? InS = true : InS = false;

    // Parameters calculation
    //   
    double A00 = A[0], A11=A[1], A22=A[2];
    R[0]+=(A2+A3+A4)*S3; A[0] = ((A0+2.*A3)+(A5+A6))*C33;
    R[1]+=(B2+B3+B4)*S3; A[1] = ((B0+2.*B3)+(B5+B6))*C33;
    R[2]+=(C2+C3+C4)*S3; A[2] = ((C0+2.*C3)+(C5+C6))*C33;
    double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    A[0]*=CBA; A[1]*=CBA; A[2]*=CBA;
 
    double Sl = 2./S;  
    sA[0] = A6*Sl; 
    sA[1] = B6*Sl;
    sA[2] = C6*Sl; 

    // Jacobian calculation
    //
    for(int i=3; i!=12; ++i) {H0[i]=g0[i-3]*PS2; H1[i]=g1[i-3]*PS2; H2[i]=g2[i-3]*PS2;}

    for(int i=7; i<35; i+=7) {

      double* dR   = &(pCache.pVector[i])                           ;  
      double* dA   = &(pCache.pVector[i+3])                         ;
      double dH0   = H0[ 3]*dR[0]+H0[ 4]*dR[1]+H0[ 5]*dR[2]         ; // dHx/dp
      double dH1   = H0[ 6]*dR[0]+H0[ 7]*dR[1]+H0[ 8]*dR[2]         ; // dHy/dp
      double dH2   = H0[ 9]*dR[0]+H0[10]*dR[1]+H0[11]*dR[2]         ; // dHz/dp
      double dA0   =(H0[ 2]*dA[1]-H0[ 1]*dA[2])+(A[1]*dH2-A[2]*dH1) ; // dA0/dp
      double dB0   =(H0[ 0]*dA[2]-H0[ 2]*dA[0])+(A[2]*dH0-A[0]*dH2) ; // dB0/dp
      double dC0   =(H0[ 1]*dA[0]-H0[ 0]*dA[1])+(A[0]*dH1-A[1]*dH0) ; // dC0/dp
      double dA2   = dA0+dA[0], dX = dR[0]+(dA2+dA[0])*S4           ; // dX /dp
      double dB2   = dB0+dA[1], dY = dR[1]+(dB2+dA[1])*S4           ; // dY /dp
      double dC2   = dC0+dA[2], dZ = dR[2]+(dC2+dA[2])*S4           ; // dZ /dp
      dH0          = H1[ 3]*dX   +H1[ 4]*dY   +H1[ 5]*dZ            ; // dHx/dp
      dH1          = H1[ 6]*dX   +H1[ 7]*dY   +H1[ 8]*dZ            ; // dHy/dp
      dH2          = H1[ 9]*dX   +H1[10]*dY   +H1[11]*dZ            ; // dHz/dp
      double dA3   =(dA[0]+dB2*H1[2]-dC2*H1[1])+(B2*dH2-C2*dH1)     ; // dA3/dp
      double dB3   =(dA[1]+dC2*H1[0]-dA2*H1[2])+(C2*dH0-A2*dH2)     ; // dB3/dp
      double dC3   =(dA[2]+dA2*H1[1]-dB2*H1[0])+(A2*dH1-B2*dH0)     ; // dC3/dp
      double dA4   =(dA[0]+dB3*H1[2]-dC3*H1[1])+(B3*dH2-C3*dH1)     ; // dA4/dp
      double dB4   =(dA[1]+dC3*H1[0]-dA3*H1[2])+(C3*dH0-A3*dH2)     ; // dB4/dp
      double dC4   =(dA[2]+dA3*H1[1]-dB3*H1[0])+(A3*dH1-B3*dH0)     ; // dC4/dp
      double dA5   = dA4+dA4-dA[0];  dX = dR[0]+dA4*S               ; // dX /dp 
      double dB5   = dB4+dB4-dA[1];  dY = dR[1]+dB4*S               ; // dY /dp
      double dC5   = dC4+dC4-dA[2];  dZ = dR[2]+dC4*S               ; // dZ /dp
      dH0          = H2[ 3]*dX   +H2[ 4]*dY   +H2[ 5]*dZ            ; // dHx/dp
      dH1          = H2[ 6]*dX   +H2[ 7]*dY   +H2[ 8]*dZ            ; // dHy/dp
      dH2          = H2[ 9]*dX   +H2[10]*dY   +H2[11]*dZ            ; // dHz/dp
      double dA6   =(dB5*H2[2]-dC5*H2[1])+(B5*dH2-C5*dH1)           ; // dA6/dp
      double dB6   =(dC5*H2[0]-dA5*H2[2])+(C5*dH0-A5*dH2)           ; // dB6/dp
      double dC6   =(dA5*H2[1]-dB5*H2[0])+(A5*dH1-B5*dH0)           ; // dC6/dp
      dR[0]+=(dA2+dA3+dA4)*S3; dA[0]=((dA0+2.*dA3)+(dA5+dA6))*C33   ;      
      dR[1]+=(dB2+dB3+dB4)*S3; dA[1]=((dB0+2.*dB3)+(dB5+dB6))*C33   ; 
      dR[2]+=(dC2+dC3+dC4)*S3; dA[2]=((dC0+2.*dC3)+(dC5+dC6))*C33   ;
    }

    double* dR   = &(pCache.pVector[35])                           ;  
    double* dA   = &(pCache.pVector[38])                           ;
      
    double dH0   = H0[ 3]*dR[0]+H0[ 4]*dR[1]+H0[ 5]*dR[2]                ; // dHx/dp
    double dH1   = H0[ 6]*dR[0]+H0[ 7]*dR[1]+H0[ 8]*dR[2]                ; // dHy/dp
    double dH2   = H0[ 9]*dR[0]+H0[10]*dR[1]+H0[11]*dR[2]                ; // dHz/dp
    double dA0   =(H0[ 2]*dA[1]-H0[ 1]*dA[2])+(A[1]*dH2-A[2]*dH1+A0)     ; // dA0/dp
    double dB0   =(H0[ 0]*dA[2]-H0[ 2]*dA[0])+(A[2]*dH0-A[0]*dH2+B0)     ; // dB0/dp
    double dC0   =(H0[ 1]*dA[0]-H0[ 0]*dA[1])+(A[0]*dH1-A[1]*dH0+C0)     ; // dC0/dp
    double dA2   = dA0+dA[0], dX = dR[0]+(dA2+dA[0])*S4                  ; // dX /dp
    double dB2   = dB0+dA[1], dY = dR[1]+(dB2+dA[1])*S4                  ; // dY /dp
    double dC2   = dC0+dA[2], dZ = dR[2]+(dC2+dA[2])*S4                  ; // dZ /dp
    dH0          = H1[ 3]*dX   +H1[ 4]*dY   +H1[ 5]*dZ                   ; // dHx/dp
    dH1          = H1[ 6]*dX   +H1[ 7]*dY   +H1[ 8]*dZ                   ; // dHy/dp
    dH2          = H1[ 9]*dX   +H1[10]*dY   +H1[11]*dZ                   ; // dHz/dp
    double dA3   =(dA[0]+dB2*H1[2]-dC2*H1[1])+((B2*dH2-C2*dH1)+(A3-A00)) ; // dA3/dp
    double dB3   =(dA[1]+dC2*H1[0]-dA2*H1[2])+((C2*dH0-A2*dH2)+(B3-A11)) ; // dB3/dp
    double dC3   =(dA[2]+dA2*H1[1]-dB2*H1[0])+((A2*dH1-B2*dH0)+(C3-A22)) ; // dC3/dp
    double dA4   =(dA[0]+dB3*H1[2]-dC3*H1[1])+((B3*dH2-C3*dH1)+(A4-A00)) ; // dA4/dp
    double dB4   =(dA[1]+dC3*H1[0]-dA3*H1[2])+((C3*dH0-A3*dH2)+(B4-A11)) ; // dB4/dp
    double dC4   =(dA[2]+dA3*H1[1]-dB3*H1[0])+((A3*dH1-B3*dH0)+(C4-A22)) ; // dC4/dp
    double dA5   = dA4+dA4-dA[0];  dX = dR[0]+dA4*S                      ; // dX /dp 
    double dB5   = dB4+dB4-dA[1];  dY = dR[1]+dB4*S                      ; // dY /dp
    double dC5   = dC4+dC4-dA[2];  dZ = dR[2]+dC4*S                      ; // dZ /dp
    dH0          = H2[ 3]*dX   +H2[ 4]*dY   +H2[ 5]*dZ                   ; // dHx/dp
    dH1          = H2[ 6]*dX   +H2[ 7]*dY   +H2[ 8]*dZ                   ; // dHy/dp
    dH2          = H2[ 9]*dX   +H2[10]*dY   +H2[11]*dZ                   ; // dHz/dp
    double dA6   =(dB5*H2[2]-dC5*H2[1])+(B5*dH2-C5*dH1+A6)               ; // dA6/dp
    double dB6   =(dC5*H2[0]-dA5*H2[2])+(C5*dH0-A5*dH2+B6)               ; // dB6/dp
    double dC6   =(dA5*H2[1]-dB5*H2[0])+(A5*dH1-B5*dH0+C6)               ; // dC6/dp
    dR[0]+=(dA2+dA3+dA4)*S3; dA[0]=((dA0+2.*dA3)+(dA5+dA6))*C33          ;      
    dR[1]+=(dB2+dB3+dB4)*S3; dA[1]=((dB0+2.*dB3)+(dB5+dB6))*C33          ; 
    dR[2]+=(dC2+dC3+dC4)*S3; dA[2]=((dC0+2.*dC3)+(dC5+dC6))*C33          ;
    return S;
  }
  return S;
}

/////////////////////////////////////////////////////////////////////////////////
// Test new cross point
/////////////////////////////////////////////////////////////////////////////////
bool Acts::RungeKuttaEngine::newCrossPoint(const CylinderSurface& Su, const double* Ro, const double* P ) const
{
  const double pi = 3.1415927, pi2=2.*pi; 
  const Transform3D& T = Su.transform();
  double Ax[3] = {T(0,0),T(1,0),T(2,0)};
  double Ay[3] = {T(0,1),T(1,1),T(2,1)};

  double R     = Su.bounds().r();
  double x     = Ro[0]-T(0,3);
  double y     = Ro[1]-T(1,3);
  double z     = Ro[2]-T(2,3);

  double RC    = x*Ax[0]+y*Ax[1]+z*Ax[2];
  double RS    = x*Ay[0]+y*Ay[1]+z*Ay[2];

  if( (RC*RC+RS*RS) <= (R*R) ) return false;
  
  x           = P[0]-T(0,3);
  y           = P[1]-T(1,3);
  z           = P[2]-T(2,3);
  RC          = x*Ax[0]+y*Ax[1]+z*Ax[2];
  RS          = x*Ay[0]+y*Ay[1]+z*Ay[2];
  double dF   = fabs(atan2(RS,RC)-Su.bounds().averagePhi());
  if(dF > pi) dF = pi2-pi;
  if(dF <= Su.bounds().halfPhiSector()) return false;
  return true;
}

/////////////////////////////////////////////////////////////////////////////////
// Straight line trajectory model 
/////////////////////////////////////////////////////////////////////////////////
double Acts::RungeKuttaEngine::straightLineStep(int navigationStep, PropagationCache& pCache, double  S) const
{
    
  EX_MSG_VERBOSE(navigationStep, "propagate", "<T> ", "straightLineStep called"); 
    
  double*  R   = &(pCache.pVector[ 0]);             // Start coordinates
  double*  A   = &(pCache.pVector[ 3]);             // Start directions
  double* sA   = &(pCache.pVector[42]);           

  // Track parameters in last point
  R[0]+=(A[0]*S); R[1]+=(A[1]*S); R[2]+=(A[2]*S); if(!pCache.useJacobian) return S;
  
  // Derivatives of track parameters in last point
  for(int i=7; i<42; i+=7) {

    double* dR = &(pCache.pVector[i]); 
    double* dA = &(pCache.pVector[i+3]);
    dR[0]+=(dA[0]*S); dR[1]+=(dA[1]*S); dR[2]+=(dA[2]*S);
  }
  sA[0]=sA[1]=sA[2]=0.; return S;
}

/////////////////////////////////////////////////////////////////////////////////
// Build new track parameters without propagation
/////////////////////////////////////////////////////////////////////////////////
const Acts::TrackParameters* Acts::RungeKuttaEngine::buildTrackParametersWithoutPropagation(const TrackParameters& tParameters,double* jacobian) const
{
  jacobian[0]=jacobian[6]=jacobian[12]=jacobian[18]=jacobian[20]=1.;
  jacobian[1]=jacobian[2]=jacobian[3]=jacobian[4]=jacobian[5]=jacobian[7]=jacobian[8]=jacobian[9]=jacobian[10]=jacobian[11]=jacobian[13]=jacobian[14]=jacobian[15]=jacobian[16]=jacobian[17]=jacobian[19]=0.;
  return tParameters.clone();
}

/////////////////////////////////////////////////////////////////////////////////
// Build new neutral track parameters without propagation
/////////////////////////////////////////////////////////////////////////////////
const Acts::NeutralParameters* Acts::RungeKuttaEngine::buildNeutralParametersWithoutPropagation(const NeutralParameters& nParameters,double* jacobian) const
{
  jacobian[0]=jacobian[6]=jacobian[12]=jacobian[18]=jacobian[20]=1.;
  jacobian[1]=jacobian[2]=jacobian[3]=jacobian[4]=jacobian[5]=jacobian[7]=jacobian[8]=jacobian[9]=jacobian[10]=jacobian[11]=jacobian[13]=jacobian[14]=jacobian[15]=jacobian[16]=jacobian[17]=jacobian[19]=0.;
  return nParameters.clone();
}

/////////////////////////////////////////////////////////////////////////////////
// Step estimator take into accout curvature
/////////////////////////////////////////////////////////////////////////////////
double Acts::RungeKuttaEngine::stepEstimatorWithCurvature(PropagationCache& pCache, int kind, double* Su, bool& Q) const
{
  // Straight step estimation
  //
  RungeKuttaUtils utils;
  
  double  Step = utils.stepEstimator(kind,Su,pCache.pVector,Q); if(!Q) return 0.; 
  double AStep = fabs(Step);
  if ( kind || AStep < m_straightStep || !pCache.mcondition ) return Step;

  const double* SA = &(pCache.pVector[42]); // Start direction
  double S = .5*Step;
  
  double Ax    = pCache.pVector[3]+S*SA[0];
  double Ay    = pCache.pVector[4]+S*SA[1];
  double Az    = pCache.pVector[5]+S*SA[2];
  double As    = 1./sqrt(Ax*Ax+Ay*Ay+Az*Az);

  double PN[6] = {pCache.pVector[0],pCache.pVector[1],pCache.pVector[2],Ax*As,Ay*As,Az*As};
  double StepN = utils.stepEstimator(kind,Su,PN,Q); if(!Q) {Q = true; return Step;}
  if(fabs(StepN) < AStep) return StepN; return Step;
} 

