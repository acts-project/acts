//////////////////////////////////////////////////////////////////
// ExtrapolationEngineTest.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Test module
#include "ExtrapolationTest/ExtrapolationEngineTest.h"
// Geometry module
#include "Surfaces/PerigeeSurface.h"
// Gaudi
#include "GaudiKernel/ITHistSvc.h"

DECLARE_COMPONENT(Acts::ExtrapolationEngineTest)

Acts::ExtrapolationEngineTest::ExtrapolationEngineTest(const std::string& name, ISvcLocator* pSvcLocator) :
 Acts::ExtrapolationTestBase(name, pSvcLocator),
 m_extrapolationEngine("",name),
 m_parametersProcessor(""),
 m_processSensitive(true),
 m_processPassive(true),
 m_parametersMode(1),
 m_particleHypothesis(2),
 m_smearProductionVertex(false),
 m_smearFlatOriginT(false),
 m_smearFlatOriginZ(false),
 m_sigmaOriginT(0.),
 m_sigmaOriginZ(0.),
 m_d0Min(0.),
 m_d0Max(0.),
 m_z0Min(0.),
 m_z0Max(0.),
 m_etaMin(-3.),
 m_etaMax(3.),
 m_phiMin(-M_PI),
 m_phiMax(M_PI),
 m_ptMin(100.),
 m_ptMax(100000.),
 m_pathLimit(10e10),
 m_collectSensitive(false),
 m_collectPassive(false),
 m_collectBoundary(false),
 m_collectMaterial(false),
 m_collectJacobians(false),
 m_sensitiveCurvilinear(false),
 m_searchMode(0), 
 m_backExtrapolation(false),
 m_stepsPhi(1),
 m_currentPhiStep(0),
 m_etaScans(0),
 m_currentEta(0.),
 m_phiScans(0),
 m_currentPhi(0.),
 m_splitCharge(false),
 m_writeTree(true),
 m_treeName("ExtrapolationEngineTest"),
 m_treeFolder("/val/"),
 m_treeDescription("ExtrapolationEngine test setup"),
 m_tree(0),
 m_tRandom(0),
 m_startPositionX(0),
 m_startPositionY(0),
 m_startPositionZ(0),
 m_startPositionR(0),
 m_startPhi(0),
 m_startTheta(0),
 m_startEta(0),
 m_startP(0),
 m_startPt(0),
 m_charge(-1.0),
 m_endSuccessful(0),
 m_endPositionX(0),
 m_endPositionY(0),
 m_endPositionZ(0),
 m_endPositionR(0),
 m_endPhi(0),
 m_endTheta(0),
 m_endEta(0),
 m_endP(0),
 m_endPt(0),
 m_endPathLength(0),
 m_backSuccessful(0),
 m_backPositionX(0),
 m_backPositionY(0),
 m_backPositionZ(0),
 m_backPositionR(0),
 m_backPhi(0),
 m_backTheta(0),
 m_backEta(0),
 m_backP(0),
 m_backPt(0),
 m_sensitiveLocalType(nullptr),
 m_sensitiveLocal0(nullptr),
 m_sensitiveLocal1(nullptr),
 m_materialThicknessInX0Accumulated(nullptr),
 m_materialThicknessInX0Steps(nullptr),
 m_materialThicknessInL0Steps(nullptr),
 m_materialPositionX(nullptr),
 m_materialPositionY(nullptr),
 m_materialPositionZ(nullptr),
 m_materialPositionR(nullptr),
 m_materialPositionP(nullptr),
 m_materialPositionPt(nullptr),
 m_materialScaling(nullptr)
{
    // the extrapolation engine
    declareProperty("ExtrapolationEngine",      m_extrapolationEngine);
    // the Parameter Base
    declareProperty("ParametersProcessor",      m_parametersProcessor);
    declareProperty("ProcessSensitive",         m_processSensitive);
    declareProperty("ProcessPassive",           m_processPassive);
    // charged / neutral & other stuff
    declareProperty("ParametersMode",           m_parametersMode);
    declareProperty("ParticleHypothesis",       m_particleHypothesis);
    declareProperty("BackExtrapolation",        m_backExtrapolation);
    // configuration
    declareProperty("PathLimit",                m_pathLimit);
    declareProperty("CollectSensitive",         m_collectSensitive);
    declareProperty("CollectPassive",           m_collectPassive);
    declareProperty("CollectBoundary",          m_collectBoundary);
    declareProperty("CollectMaterial",          m_collectMaterial);
    declareProperty("CollectJacobians",         m_collectJacobians);
    declareProperty("SensitiveCurvilinear",     m_sensitiveCurvilinear);
    declareProperty("SearchMode",               m_searchMode);
    // Mode for scanning in steps
    declareProperty("ScanMode",                 m_scanMode);
    declareProperty("EtaScans",                 m_etaScans);
    declareProperty("PhiSteps",                 m_stepsPhi);
    declareProperty("PhiScans",                 m_phiScans);
    declareProperty("SplitCharge",              m_splitCharge);
    // production vertices
    declareProperty("SmearOrigin",              m_smearProductionVertex);
    declareProperty("SmearFlatOriginD0",        m_smearFlatOriginT);
    declareProperty("SmearFlatOriginZ0",        m_smearFlatOriginZ);
    declareProperty("SimgaOriginD0",            m_sigmaOriginT);
    declareProperty("SimgaOriginZ0",            m_sigmaOriginZ);
    // d0 min / max values for flat smearing
    declareProperty("D0Min",                    m_d0Min);
    declareProperty("D0Max",                    m_d0Max);
    // z0 min / max values for flat smearing
    declareProperty("Z0Min",                    m_z0Min);
    declareProperty("Z0Max",                    m_z0Max);
    // eta min / max values
    declareProperty("EtaMin",                   m_etaMin);
    declareProperty("EtaMax",                   m_etaMax);
    // phi min / max values
    declareProperty("PhiMin",                   m_phiMin);
    declareProperty("PhiMax",                   m_phiMax);
    // pt min / max values
    declareProperty("PtMin",                    m_ptMin);
    declareProperty("PtMax",                    m_ptMax);
    // the properties
    declareProperty("NumTestsPerEvent",         m_numTests=100);
    declareProperty("WriteValidationTree",      m_writeTree);    
    declareProperty("TreeName",                 m_treeName);
    declareProperty("TreeFolder",               m_treeFolder);
    declareProperty("TreeDescription",          m_treeDescription);
}

StatusCode Acts::ExtrapolationEngineTest::finalize() {

    if (m_writeTree){
        // memory clean up
        for (size_t ip = 0; ip < m_parameterNames.size(); ++ip){
            // create
            delete m_pPositionX[ip];
            delete m_pPositionY[ip];
            delete m_pPositionZ[ip];
            delete m_pPositionR[ip];
            delete m_pPhi[ip];
            delete m_pTheta[ip];
            delete m_pEta[ip];
            delete m_pP[ip];
            delete m_pPt[ip];
        }  
        if (m_collectSensitive){
            delete m_sensitiveLocalType;
            delete m_sensitiveLocal0;
            delete m_sensitiveLocal1;
        }
        if (m_collectMaterial){
            delete m_materialThicknessInX0Accumulated;
            delete m_materialThicknessInX0Steps;
            delete m_materialThicknessInL0Steps;
            delete m_materialPositionX;
            delete m_materialPositionY;
            delete m_materialPositionZ;
            delete m_materialPositionR;
            delete m_materialPositionP;
            delete m_materialPositionPt;
            delete m_materialScaling;
        }
    }
    return StatusCode::SUCCESS;
}

StatusCode Acts::ExtrapolationEngineTest::initializeTest()
{
    MSG_INFO("initialize ExtrapolationEngineTest0");
    // Extrapolation engine
    RETRIEVE_FATAL(m_extrapolationEngine);
    MSG_INFO("initialize ExtrapolationEngineTest1");
    // Writer
    RETRIEVE_NONEMPTY_FATAL(m_parametersProcessor);
    // success
    return StatusCode::SUCCESS;
}

StatusCode Acts::ExtrapolationEngineTest::bookTree()
{
    MSG_VERBOSE("Booking the Extrapolation test Tree.");
    
    if (m_writeTree){
        // ------------------------------> OUTPUT NTUPLE (geometry validation)
        m_tree = new TTree(m_treeName.c_str(), m_treeDescription.c_str());
        // add the Branches
        m_tree->Branch("StartPosX",     &m_startPositionX);
        m_tree->Branch("StartPosY",     &m_startPositionY);
        m_tree->Branch("StartPosZ",     &m_startPositionZ);
        m_tree->Branch("StartPosR",     &m_startPositionR);
        m_tree->Branch("StartPhi",      &m_startPhi);
        m_tree->Branch("StartEta",      &m_startEta);
        m_tree->Branch("StartTheta",    &m_startTheta);
        m_tree->Branch("StartP",        &m_startP);
        m_tree->Branch("StartPt",       &m_startPt);
        if (m_parametersMode) m_tree->Branch("StartCharge",   &m_charge);
        
        m_tree->Branch("EndSuccessful", &m_endSuccessful);
        m_tree->Branch("EndPosX",       &m_endPositionX);
        m_tree->Branch("EndPosY",       &m_endPositionY);
        m_tree->Branch("EndPosZ",       &m_endPositionZ);
        m_tree->Branch("EndPosR",       &m_endPositionR);
        m_tree->Branch("EndPhi",        &m_endPhi);
        m_tree->Branch("EndEta",        &m_endEta);
        m_tree->Branch("EndTheta",      &m_endTheta);
        m_tree->Branch("EndP",          &m_endP);
        m_tree->Branch("EndPt",         &m_endPt);
        // also add the path length
        m_tree->Branch("EndPathLength", &m_endPathLength);
        
        if (m_backExtrapolation){
            m_tree->Branch("BackSuccessful", &m_backSuccessful);
            m_tree->Branch("BackPosX",       &m_backPositionX);
            m_tree->Branch("BackPosY",       &m_backPositionY);
            m_tree->Branch("BackPosZ",       &m_backPositionZ);
            m_tree->Branch("BackPosR",       &m_backPositionR);
            m_tree->Branch("BackPhi",        &m_backPhi);
            m_tree->Branch("BackEta",        &m_backEta);
            m_tree->Branch("BackTheta",      &m_backTheta);
            m_tree->Branch("BackP",          &m_backP);
            m_tree->Branch("BackPt",         &m_backPt);
        }
        
        // this fixes the parameters to order
        m_parameterNames.push_back("Sensitive");
        m_parameterNames.push_back("Passive");
        m_parameterNames.push_back("Boundary");
        for (size_t ip = 0; ip < m_parameterNames.size(); ++ip){
            // create
            m_pPositionX.push_back( new std::vector<float> );
            m_pPositionY.push_back( new std::vector<float> );
            m_pPositionZ.push_back( new std::vector<float> );
            m_pPositionR.push_back( new std::vector<float> );
            m_pPhi.push_back( new std::vector<float> );      
            m_pTheta.push_back( new std::vector<float> );
            m_pEta.push_back( new std::vector<float> );
            m_pP.push_back( new std::vector<float> );
            m_pPt.push_back( new std::vector<float> );
	    
            // define the branches    
            m_tree->Branch(m_parameterNames[ip]+"PosX",       m_pPositionX[ip]);
            m_tree->Branch(m_parameterNames[ip]+"PosY",       m_pPositionY[ip]);
            m_tree->Branch(m_parameterNames[ip]+"PosZ",       m_pPositionZ[ip]);
            m_tree->Branch(m_parameterNames[ip]+"PosR",       m_pPositionR[ip]);
            m_tree->Branch(m_parameterNames[ip]+"Phi",        m_pPhi[ip]      );
            m_tree->Branch(m_parameterNames[ip]+"Eta",        m_pTheta[ip]    );
            m_tree->Branch(m_parameterNames[ip]+"Theta",      m_pEta[ip]      );
            m_tree->Branch(m_parameterNames[ip]+"P",          m_pP[ip]        );
            m_tree->Branch(m_parameterNames[ip]+"Pt",         m_pPt[ip]       );
        }
        
        if (m_collectSensitive){
             m_sensitiveLocalType              = new std::vector<float>;
             m_sensitiveLocal0                 = new std::vector<float>;
             m_sensitiveLocal1                 = new std::vector<float>;
             m_tree->Branch("SensitiveType",   m_sensitiveLocalType);
             m_tree->Branch("SensitiveLoc0",   m_sensitiveLocal0);
             m_tree->Branch("SensitiveLoc1",   m_sensitiveLocal1);
        }
        
        // collect the material, you need branches for this
        if (m_collectMaterial){
            m_materialThicknessInX0Accumulated = new std::vector<float>;
            m_materialThicknessInX0Steps       = new std::vector<float>;
            m_materialThicknessInL0Steps       = new std::vector<float>;     
            m_materialPositionX                = new std::vector<float>;
            m_materialPositionY                = new std::vector<float>;
            m_materialPositionZ                = new std::vector<float>;
            m_materialPositionR                = new std::vector<float>;
            m_materialPositionP                = new std::vector<float>;
            m_materialPositionPt               = new std::vector<float>;
            m_materialScaling                  = new std::vector<float>;   
            m_tree->Branch("MaterialThicknessInX0",                &m_materialThicknessInX0);
            m_tree->Branch("MaterialThicknessInL0",                &m_materialThicknessInL0);
            m_tree->Branch("MaterialThicknessZARho",               &m_materialThicknessZARho);
            m_tree->Branch("MaterialThicknessSensitiveInX0",       &m_materialThicknessInX0Sensitive);
            m_tree->Branch("MaterialThicknessPassiveInX0",         &m_materialThicknessInX0Passive  );
            m_tree->Branch("MaterialThicknessBoundaryInX0",        &m_materialThicknessInX0Boundary );
            m_tree->Branch("MaterialThicknessAccumulatedX0",       m_materialThicknessInX0Accumulated );
            m_tree->Branch("MaterialThicknessStepsInX0",           m_materialThicknessInX0Steps);
            m_tree->Branch("MaterialThicknessStepsInL0",           m_materialThicknessInL0Steps);
            m_tree->Branch("MaterialPosX",                         m_materialPositionX);
            m_tree->Branch("MaterialPosY",                         m_materialPositionY);
            m_tree->Branch("MaterialPosZ",                         m_materialPositionZ);
            m_tree->Branch("MaterialPosR",                         m_materialPositionR);
            m_tree->Branch("MaterialPosP",                         m_materialPositionP);
            m_tree->Branch("MaterialPosPt",                        m_materialPositionPt);       
            m_tree->Branch("MaterialScaling",                      m_materialScaling);
        }
        
        // now register the Tree
        ITHistSvc* tHistSvc = 0;
        if (service("THistSvc",tHistSvc).isFailure()) {
          MSG_ERROR( "initialize() Could not find Hist Service  -> Switching Tree output off !" );
            delete m_tree; m_tree = 0;
        }
        if (tHistSvc && ((tHistSvc->regTree(m_treeFolder+m_treeName, m_tree)).isFailure()) ) {
            MSG_ERROR( "initialize() Could not register the validation Tree -> Switching Tree output off !" );
            delete m_tree; m_tree = 0;
        }   
    } 
    return StatusCode::SUCCESS;

}

StatusCode Acts::ExtrapolationEngineTest::runTest()
{
  MSG_VERBOSE("Running the ExtrapolationEngineTest Test in parameters mode : " << m_parametersMode);

  if (!m_parametersProcessor.empty() &&  m_parametersProcessor->initProcessor().isFailure() )
      MSG_WARNING("Problem initializing the Processor");


  // ----------------- creation of the surfaces & track parameters -------------
  for (size_t it = 0; it < ExtrapolationTestBase::m_numTests; ++it ){
      // verbose output
      MSG_DEBUG("===> starting test " << it << " <<===");
      // create the curvilinear parameters
      double eta   = (m_scanMode) ? m_currentEta : m_etaMin + (m_etaMax-m_etaMin)*Acts::ExtrapolationTestBase::m_flatDist->shoot();
      double theta = 2.*atan(exp(-eta));
      double phi   = (m_scanMode) ? m_currentPhi : m_phiMin + (m_phiMax-m_phiMin)*Acts::ExtrapolationTestBase::m_flatDist->shoot();
      double pt    = m_ptMin  + (m_ptMax - m_ptMin)*Acts::ExtrapolationTestBase::m_flatDist->shoot();
      double p     = pt/sin(theta);
      double q     = m_splitCharge ? m_charge*-1. : (m_parametersMode ? (Acts::ExtrapolationTestBase::m_flatDist->shoot() > 0.5 ? 1. : -1) : 1.);      // charge or neutral

      // initializa the validation variables
      m_endSuccessful= 0;
      m_endPositionX = 0.;
      m_endPositionY = 0.;
      m_endPositionZ = 0.;
      m_endPositionR = 0.;
      m_endPhi       = 0.;
      m_endEta       = 0.;
      m_endTheta     = 0.;
      m_endP         = 0.;
      m_endPt        = 0.;
      
      if (m_writeTree){
        for (size_t ip = 0; ip < m_parameterNames.size(); ++ip){
            // clear
            m_pPositionX[ip]->clear();
            m_pPositionY[ip]->clear();
            m_pPositionZ[ip]->clear();
            m_pPositionR[ip]->clear();
            m_pPhi[ip]->clear();    
            m_pTheta[ip]->clear();
            m_pEta[ip]->clear();
            m_pP[ip]->clear();
            m_pPt[ip]->clear();            
        }  
        
        if (m_collectSensitive){
            m_sensitiveLocalType->clear();
            m_sensitiveLocal0->clear();
            m_sensitiveLocal1->clear();
        }
        
        // material collection
        m_materialThicknessInX0                 = 0.;
        m_materialThicknessInL0                 = 0.;
        m_materialThicknessZARho                = 0.;
        m_materialThicknessInX0Sensitive        = 0.;
        m_materialThicknessInX0Passive          = 0.;
        m_materialThicknessInX0Boundary         = 0.;
        if (m_collectMaterial){
            m_materialThicknessInX0Accumulated->clear();
            m_materialThicknessInX0Steps->clear();
            m_materialThicknessInX0Steps->clear();
            m_materialPositionX->clear();
            m_materialPositionY->clear();
            m_materialPositionZ->clear();
            m_materialPositionR->clear();
            m_materialPositionP->clear();
            m_materialPositionPt->clear();
            m_materialScaling->clear();
        }
      }              
          
      Vector3D momentum(p*sin(theta)*cos(phi), p*sin(theta)*sin(phi), p*cos(theta));        

      // create the start parameters
      double d0 = m_smearProductionVertex ? (m_smearFlatOriginT ? (m_d0Min + (m_d0Max-m_d0Min)*Acts::ExtrapolationTestBase::m_flatDist->shoot()) : Acts::ExtrapolationTestBase::m_gaussDist->shoot()*m_sigmaOriginT) : 0.;
      double z0 = m_smearProductionVertex ? (m_smearFlatOriginZ ? (m_z0Min + (m_z0Max-m_z0Min)*Acts::ExtrapolationTestBase::m_flatDist->shoot()) : Acts::ExtrapolationTestBase::m_gaussDist->shoot()*m_sigmaOriginZ) : 0.;

      m_startPhi       = phi;
      m_startEta       = eta;
      m_startTheta     = theta;
      m_startPt        = pt;
      m_startP         = p;
      m_charge         = q;

            // preps
      std::unique_ptr<ActsSymMatrixD<5> > cov;
      ActsVectorD<5> pars; pars << d0, z0, phi, theta, q/p;

      MSG_VERBOSE("Building parameters from Perigee with (" << d0 << ", " << z0 << ", " << phi << ", " << theta << ", " << q/p);

      Acts::PerigeeSurface pSurface(Vector3D(0.,0.,0.));
      if (m_parametersMode == 0 ){
          // create the neutral parameters
          NeutralBoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
          // Screen output
          if (executeTestT<Acts::NeutralParameters>(startParameters).isFailure())
              MSG_WARNING("Test with neutral parameters did not succeed.");
          else if (m_writeTree)
              m_tree->Fill();
      } else {
          // create the charged parameters
          BoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
          if (executeTestT<Acts::TrackParameters>(startParameters).isFailure())
              MSG_WARNING("Test with neutral parameters did not succeed.");
          else if (m_writeTree)
              m_tree->Fill();
      }

  } // loop over tests
  return StatusCode::SUCCESS;
}

