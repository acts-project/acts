///////////////////////////////////////////////////////////////////
// PropagationEngine.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// STL
#include <sstream>
// Extrapolation module
#include "ExtrapolationEngine/PropagationEngine.h"
#include "ExtrapolationInterfaces/IPropagator.h"
// Geometry module
#include "Surfaces/Surface.h"
// EventData module
#include "EventDataUtils/PropDirection.h"


DECLARE_SERVICE_FACTORY(Acts::PropagationEngine)

// constructor
Acts::PropagationEngine::PropagationEngine(const std::string& name, ISvcLocator* svc):
  Acts::ServiceBase(name, svc),
  m_propagator(""),
  m_pathLimitTolerance(0.01)
{
    // configure the Propagator
    declareProperty("Propagator"             , m_propagator);
    // steering of the screen outoput (SOP)
    declareProperty("OutputPrefix"           , m_sopPrefix);
    declareProperty("OutputPostfix"          , m_sopPostfix);
    // the path limit tolerance
    declareProperty("PathLimitTolerance"     , m_pathLimitTolerance);
}

// destructor
Acts::PropagationEngine::~PropagationEngine()
{}

/** Query the interfaces. */
StatusCode Acts::PropagationEngine::queryInterface(const InterfaceID& riid, void** ppvInterface)
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

// the interface method initialize
StatusCode Acts::PropagationEngine::initialize()
{
    MSG_DEBUG("initialize()" );
    //Service needs to be initialized
    if (!ServiceBase::initialize()) return StatusCode::FAILURE;
    RETRIEVE_NONEMPTY_FATAL(m_propagator);
    return StatusCode::SUCCESS;
}

// the interface method finalize
StatusCode Acts::PropagationEngine::finalize()
{    
    MSG_DEBUG("finalize()" );
    return StatusCode::SUCCESS;
}

/** propagation method for charged particles */
Acts::ExtrapolationCode Acts::PropagationEngine::propagate(Acts::ExCellCharged& eCell,
                                                         const Acts::Surface& sf,
                                                         Acts::PropDirection pDir,
                                                         const Acts::BoundaryCheck& bcheck,
                                                         bool returnCurvilinear) const
{
    EX_MSG_DEBUG(++eCell.navigationStep, "propagate", "char", "propagation engine called with charged parameters with propagation direction " << pDir );

    double propLength = -1.;
    if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::StopWithPathLimit)){
        // the path limit
        propLength = eCell.pathLimit > 0 ? (eCell.pathLimit-eCell.pathLength) : eCell.pathLimit;
        EX_MSG_VERBOSE(eCell.navigationStep, "propagate", "char", "available step length for this propagation " << propLength );
    }
    // it is the final propagation if it is the endSurface
    bool finalPropagation = (eCell.endSurface == (&sf));

    Acts::TransportJacobian* tjac = 0;
   
    // propagate using the IPropagator
    const Acts::TrackParameters* pParameters = m_propagator->propagate(*eCell.leadParameters, 
                                                                      sf,
                                                                      pDir,
                                                                      bcheck,
                                                                      eCell.mFieldMode,
                                                                      tjac,
                                                                      propLength,
                                                                      returnCurvilinear);

    // set the return type according to how the propagation went
    if (pParameters){
       // cache the last lead parameters, useful in case a navigation error occured
       eCell.lastLeadParameters = eCell.leadParameters;
       // assign the lead and end parameters
       eCell.leadParameters = pParameters;
       // check what to do with the path Length
       if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::StopWithPathLimit) || propLength > 0){
           // add the new propagation length to the path length
           eCell.pathLength += propLength;
           // check if Limit reached
           if (eCell.pathLimitReached(m_pathLimitTolerance)){
               EX_MSG_VERBOSE(eCell.navigationStep, "propagate", "char", "path limit of " << eCell.pathLimit << " successfully reached -> stopping." ); 
               return Acts::ExtrapolationCode::SuccessPathLimit;
           }
       }
       // check if the propagation was called with directly, then lead parameters become end parameters
       if (eCell.checkConfigurationMode(Acts::ExtrapolationMode::Direct)) 
           eCell.endParameters = eCell.leadParameters;
	 
       // return Success only if it is the final propagation - the extrapolation engine knows that 
       return (finalPropagation ? Acts::ExtrapolationCode::SuccessDestination : Acts::ExtrapolationCode::InProgress);
    }                                                                      
    // return - recovered means that the leadParameters are the input ones 
    return (finalPropagation ? Acts::ExtrapolationCode::FailureDestination : Acts::ExtrapolationCode::Recovered) ;
}

/** propagation method for neutral particles */
Acts::ExtrapolationCode Acts::PropagationEngine::propagate(Acts::ExCellNeutral& eCell,
                                                         const Acts::Surface& sf,
                                                         Acts::PropDirection pDir,
                                                         const Acts::BoundaryCheck& bcheck,
                                                         bool returnCurvilinear) const
{
    EX_MSG_DEBUG(++eCell.navigationStep, "propagate", "neut", "propagation engine called with neutral parameters with propagation direction " << pDir );
    // leave this for the moment, can re replaced by an appropriate propagator call later
    if (eCell.leadParameters->covariance()){
        EX_MSG_VERBOSE(eCell.navigationStep,"propagate", "neut", "propagation of neutral parameters with covariances requested. This is not supported for the moment.");
    }
    // the pathLimit cache so far
    double cPath = eCell.pathLength;
    // it is the final propagation if it is the endSurface
    bool finalPropagation = (eCell.endSurface == (&sf));
    // intersect the surface
    Acts::Intersection sfIntersection = (pDir!=Acts::anyDirection) ? sf.intersectionEstimate(eCell.leadParameters->position(),
                                                                                               pDir*eCell.leadParameters->momentum().unit(),
                                                                                               true,
                                                                                               bcheck) :
                                                                   sf.intersectionEstimate(eCell.leadParameters->position(),
                                                                                               eCell.leadParameters->momentum().unit(),
                                                                                               false,
                                                                                               bcheck);
    // we have a valid intersection
    if (sfIntersection.valid){
        // fill the transport information - only if the propation direction is not 0 ('anyDirection')
        if (pDir!=Acts::anyDirection){
           double pLength = (sfIntersection.position-eCell.leadParameters->position()).mag(); 
           EX_MSG_VERBOSE(eCell.navigationStep,"propagate", "neut", "path length of " << pLength << " added to the extrapolation cell (limit = " << eCell.pathLimit << ")" );    
           eCell.stepTransport(sf,pLength);
        }
        // now check if it is valud it's further away than the pathLimit
        if (eCell.pathLimitReached(m_pathLimitTolerance)){
            // cache the last lead parameters
            eCell.lastLeadParameters = eCell.leadParameters;
            // create new neutral curvilinear parameters at the path limit reached
            double pDiff = eCell.pathLimit - cPath;
            Vector3D position = eCell.leadParameters->position()+pDiff*eCell.leadParameters->momentum().unit();
            eCell.leadParameters = new Acts::NeutralCurvilinearParameters(nullptr,position,eCell.leadParameters->momentum());
            EX_MSG_VERBOSE(eCell.navigationStep,"propagate", "neut", "path limit of " << eCell.pathLimit << " reached. Stopping extrapolation.");
            return Acts::ExtrapolationCode::SuccessPathLimit;
        }
		// cache the last lead parameters
		eCell.lastLeadParameters = eCell.leadParameters;
        // now exchange the lead parameters
        // create the new curvilinear paramters at the surface intersection -> if so, trigger the success
        eCell.leadParameters = returnCurvilinear ? static_cast<Acts::NeutralParameters*>(new Acts::NeutralCurvilinearParameters(nullptr,sfIntersection.position, eCell.leadParameters->momentum())):
            static_cast<Acts::NeutralParameters*>(new Acts::NeutralBoundParameters(nullptr,sfIntersection.position, eCell.leadParameters->momentum(),sf));

        // check if the propagation was called with directly, then lead parameters become end parameters
        if (eCell.checkConfigurationMode(ExtrapolationMode::Direct))
	          eCell.endParameters = eCell.leadParameters;

	      // return success for the final destination or in progress
        return (finalPropagation ? Acts::ExtrapolationCode::SuccessDestination : Acts::ExtrapolationCode::InProgress);

    } else {
        // give some screen output
        EX_MSG_VERBOSE(eCell.navigationStep,"propagate", "neut", "intersection with the surface did not succeed.");
    }                                                                     
   // return - recovered means that the leadParameters are the input ones 
   return (finalPropagation ? Acts::ExtrapolationCode::FailureDestination : Acts::ExtrapolationCode::Recovered) ;
}
