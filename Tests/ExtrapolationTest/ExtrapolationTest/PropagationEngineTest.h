//////////////////////////////////////////////////////////////////
// PropagationEngineTest.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ATAS_EXTRAPOLATIONTEST_PROPAGATIONENGINETEST_H
#define ATAS_EXTRAPOLATIONTEST_PROPAGATIONENGINETEST_H 1

// Athena & Gaudi includes
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
// Extrapolation module
#include "ExtrapolationInterfaces/IPropagationEngine.h"
// Test modules
#include "ExtrapolationTest/ExtrapolationTestBase.h"
#include "TestInterfaces/IParametersBaseProcessor.h"
// ROOT
#include "TString.h"
#include "TRandom.h"

class TTree;

namespace Acts {
                 
    /** @class PropagationEngineTest
               
        Test Algorithm to run test extrapolations with the new IExtrapolationEngine           
               
        @author Andreas.Salzburger@cern.ch, Noemi.Calace@cern.ch       
      */
      
   class PropagationEngineTest : public ExtrapolationTestBase  {
     
     public:

       /** Standard Athena-Algorithm Constructor */
       PropagationEngineTest(const std::string& name, ISvcLocator* pSvcLocator);
        
       /* finalize */
       StatusCode finalize();

       /* specify the test here */
       StatusCode runTest();

       /** initialize the test, i.e. retrieve the TrackingGeometry Svc */
       StatusCode initializeTest();
       
       /* book the TTree branches */
       StatusCode bookTree();
       
     private:
         
       template <class T> const T* executeTestT(const T& startParameters, const Surface& sf); 
                  
       /** retrieve it */
       ServiceHandle<IPropagationEngine>            m_propagationEngine;  
       
       /** output writer */
       ToolHandle<IParametersBaseProcessor>         m_parametersProcessor;   
              
       bool                                         m_parametersMode; // 0 - neutral, 1 - charged, 2 - multi
     
       std::vector<double>                          m_destinationRadii;
       std::vector<Surface*>                        m_destinationSurfaces;
       bool                                         m_emulatePlaneSurfaces;
     
       double                                       m_etaMin;
       double                                       m_etaMax;
       double                                       m_phiMin;
       double                                       m_phiMax;        
       double                                       m_ptMin;
       double                                       m_ptMax;       
       double                                       m_pathLimit;
              
       bool                                         m_backPropagation;
       bool                                         m_returnCurvilinear;
       
       //!< the tree              
       std::string                                  m_treeName;
       std::string                                  m_treeFolder;  
       std::string                                  m_treeDescription;
       TTree*                                       m_tree;
       TRandom                                      m_tRandom;        
                                                    
                                                    
       float                                        m_startPositionX;
       float                                        m_startPositionY;
       float                                        m_startPositionZ;
       float                                        m_startPositionR;
       float                                        m_startPhi;
       float                                        m_startTheta;                                                    
       float                                        m_startEta;                                                    
       float                                        m_startP;                                                    
       float                                        m_startPt;
       float                                        m_charge;                                                    
        
       int                                          m_endSuccessful;
       float                                        m_endPositionX;
       float                                        m_endPositionY;
       float                                        m_endPositionZ;
       float                                        m_endPositionR;
       float                                        m_endPhi;
       float                                        m_endTheta;                                                    
       float                                        m_endEta;                                                    
       float                                        m_endP;                                                    
       float                                        m_endPt;        
       float                                        m_endPathLength;
       
       int                                          m_backSuccessful;
       float                                        m_backPositionX;
       float                                        m_backPositionY;
       float                                        m_backPositionZ;
       float                                        m_backPositionR;
       float                                        m_backPhi;
       float                                        m_backTheta;                                                    
       float                                        m_backEta;                                                    
       float                                        m_backP;                                                    
       float                                        m_backPt; 
       

   };
} // end of namespace

// include the templated function
#include "PropagationEngineTest.icc"

#endif
