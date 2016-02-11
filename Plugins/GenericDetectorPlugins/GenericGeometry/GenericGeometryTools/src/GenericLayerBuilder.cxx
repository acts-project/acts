///////////////////////////////////////////////////////////////////
// GenericLayerBuilder.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Core module
#include "Algebra/AlgebraDefinitions.h"
#include "Algebra/AlgebraHelper.h"
#include "CoreInterfaces/MsgMacros.h"
// Geometry module
#include "GeometryUtils/BinUtility.h"
#include "GeometryUtils/BinnedArray2D.h"
#include "GeometryUtils/BinnedArray1D.h"
#include "GeometryUtils/BinnedArrayArray.h"
#include "Detector/CylinderLayer.h"
#include "Detector/DiscLayer.h"
#include "Surfaces/RadialBounds.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/PlanarBounds.h"
#include "Surfaces/RectangleBounds.h"
#include "Surfaces/TrapezoidBounds.h"
// A generic detector
#include "GenericGeometryTools/GenericLayerBuilder.h"
#include "GenericDetectorElement/DetectorElement.h"

DECLARE_COMPONENT(Agd::GenericLayerBuilder)

// constructor
Agd::GenericLayerBuilder::GenericLayerBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  Ats::AlgToolBase(t,n,p),
  m_layerIdentification(n),
  m_nLayers(nullptr),
  m_cLayers(nullptr),    
  m_pLayers(nullptr)   
{
    declareInterface<ILayerBuilder>(this);
   
    //  layer identificaiton
    declareProperty("LayerIdentification",          m_layerIdentification);
   
    // the central layers 
    declareProperty("CentralLayerRadii",            m_centralLayerRadii);
    declareProperty("CentralLayerEnvelopeZ",        m_centralLayerEnvelopeZ);
    declareProperty("CentralLayerModulesPhi",       m_centralModulesPhi);
    declareProperty("CentralLayerMoudlesTiltPhi",   m_centralModulesTiltPhi);        
    declareProperty("CentralLayerModulesPositionZ", m_centralModulesPositionZ);
    declareProperty("CentralLayerModuleStaggerZ",   m_centralModulesStaggerZ);
    declareProperty("CentralLayerModuleHalfX",      m_centralModuleHalfX);
    declareProperty("CentralLayerModuleHalfY",      m_centralModuleHalfY);
    declareProperty("CentralLayerModuleThickness",  m_centralModuleThickness);
    
    // the layers at p/e side 
    declareProperty("PosNegLayerPositionZ",         m_posnegLayerPositionsZ);
    declareProperty("PosNegLayerEnvelopeR",         m_posnegLayerEnvelopeR);
    declareProperty("PosNegLayerModuleRadii",       m_posnegModulesRadii);
    declareProperty("PosNegLayerModuleStaggerR",    m_posnegModuleStaggerR);        
    declareProperty("PosNegLayerModulesPhi",        m_posnegModulesPhi);
    declareProperty("PosNegLayerModulesStaggerPhi", m_posnegMoudleStaggerPhi);
    declareProperty("PosNegLayerModuleMinHalfX",    m_posnegModuleMinHalfX);
    declareProperty("PosNegLayerModuleMaxHalfX",    m_posnegModuleMaxHalfX);
    declareProperty("PosNegLayerModuleHalfY",       m_posnegModuleHalfY);
    declareProperty("PosNegLayerModuleThickness",   m_posnegModuleThickness);
    
}

// destructor
Agd::GenericLayerBuilder::~GenericLayerBuilder()
{}

// initialize
StatusCode Agd::GenericLayerBuilder::initialize()
{
    MSG_DEBUG( "initialize()" );
    return constructLayers();
}

//finalize
StatusCode Agd::GenericLayerBuilder::finalize()
{
    MSG_DEBUG( "finalize()" );
    return StatusCode::SUCCESS;
}


StatusCode Agd::GenericLayerBuilder::constructLayers() 
{
    
    // -------------------------------- central layers -----------------------------------------------------------
    typedef std::pair<const Ats::Surface*, Ats::Vector3D> SurfacePosition;
    // the central layers
    size_t numcLayers = m_centralLayerRadii.size();
    if (numcLayers){
        MSG_DEBUG("Configured to build " << numcLayers << " passive central layers.");
        m_cLayers = new Ats::LayerVector;
        m_cLayers->reserve(numcLayers);
        // loop through
        for (size_t icl = 0; icl < numcLayers; ++icl){
            // layer R/Z
            double layerR = m_centralLayerRadii[icl];
            double halfZ  = 0.;
            double minPhi = 10.;
            double maxPhi = -10.;
            // some screen output
            MSG_VERBOSE("- build layer " << icl << " with radius = " << layerR);
            // create the modules & surface array 
            Ats::SurfaceArray* sArray = nullptr;
            // surface vector 
            std::vector<SurfacePosition> sVector;
            // z/phi values for this layer
            std::vector<double> zValues   = m_centralModulesPositionZ[icl];
            std::vector<double> phiValues = m_centralModulesPhi[icl];
            std::sort(phiValues.begin(),phiValues.end());
            std::sort(zValues.begin(),zValues.end());
            // envelope stagger & cover
            double layerModleStaggerZ  = m_centralModulesStaggerZ.size() ? m_centralModulesStaggerZ[icl] : 0.;
            double layerEnvelopeCoverZ = m_centralLayerEnvelopeZ.size() ? m_centralLayerEnvelopeZ[icl] : 0.;
            double layerMinR = 10e10;
            double layerMaxR = 0;
            // module size & tilt
            double modulePhiTilt   = m_centralModulesTiltPhi[icl]; 
            double moduleHalfX     = m_centralModuleHalfX[icl];
            double moduleHalfY     = m_centralModuleHalfY[icl];
            double moduleThickness = m_centralModuleThickness[icl];
            // create the shared module 
            std::shared_ptr<const Ats::PlanarBounds> moduleBounds(new Ats::RectangleBounds(moduleHalfX,moduleHalfY));
            // now create the modules and surfaces 
            bool stagger = false;
            // Identifier @TODO unique Identifier
            size_t imodule = 0;
            sVector.reserve(zValues.size()*phiValues.size());
            // loop over z module and phi module position
            for (auto& moduleZ : zValues){
                // create the half length in z
                halfZ = halfZ > moduleZ+moduleHalfY ? halfZ :  moduleZ+moduleHalfY;
                // loop of phi values
                for (auto& modulePhi : phiValues){
                    // min/max phi
                    takeSmallerBigger(minPhi, maxPhi, modulePhi);
                    // count the modules
                    ++imodule;
                    // stagger the modules
                    double mouldeR = layerR;
                    mouldeR += stagger ? -0.5*layerModleStaggerZ : 0.5 * layerModleStaggerZ; stagger = !stagger;
                    // the position of the module
                    Ats::Vector3D moduleCenter(mouldeR*cos(modulePhi), mouldeR*sin(modulePhi), moduleZ);                     
                    // normal vectorof the surface
                    Ats::Vector3D moduleLocalZ(cos(modulePhi+modulePhiTilt),sin(modulePhi+modulePhiTilt), 0.);
                    Ats::Vector3D moduleLocalY(0.,0.,1);
                    Ats::Vector3D moduleLocalX(-sin(modulePhi+modulePhiTilt),cos(modulePhi+modulePhiTilt),0.);
                    // create the RotationMatrix
                    Ats::RotationMatrix3D moduleRotation;
                    moduleRotation.col(0) = moduleLocalX;
                    moduleRotation.col(1) = moduleLocalY;
                    moduleRotation.col(2) = moduleLocalZ;
                    // edge points 
                    double moduleEdge1 = Ats::Vector3D(moduleCenter+moduleHalfX*moduleLocalX).perp();
                    double moduleEdge2 = Ats::Vector3D(moduleCenter-moduleHalfX*moduleLocalX).perp();
                    // layer min / max 
                    takeSmallerBigger(layerMinR,layerMaxR, moduleEdge1);
                    takeSmallerBigger(layerMinR,layerMaxR, moduleEdge2);
                    // get the moduleTransform
                    std::shared_ptr<Ats::Transform3D> moduleTransform(new Ats::Transform3D(Ats::getTransformFromRotTransl(moduleRotation,moduleCenter)));
                    // create the generic detector element @TODO identifier service
                    Identifier moduleIdentifier(imodule);
                    // create the module 
                    DetectorElement* module = new DetectorElement(moduleIdentifier, moduleTransform, moduleBounds, moduleThickness);
                    // create the surface 
                    sVector.push_back(SurfacePosition(&module->surface(),moduleCenter));
                    // memory management - we need a detector store to hold them somewhere @TODO detector store facility
                    m_centralModules.push_back(module);
                }
            }
            // harmonize the phi boundaries 
            double phiStep = (maxPhi-minPhi)/(phiValues.size()-1);
            minPhi -= 0.5*phiStep;
            maxPhi += 0.5*phiStep;
            // layer thickness
            double layerThickness = (layerMaxR-layerMinR);
            // create the binUtility
            Ats::BinUtility* moduleBinUtility = new Ats::BinUtility(phiValues.size(), minPhi, maxPhi, Ats::closed, Ats::binPhi);
            (*moduleBinUtility) += Ats::BinUtility(zValues.size(), -halfZ, halfZ, Ats::open, Ats::binZ);
            // create the surface array 
            sArray = new Ats::BinnedArray2D< const Ats::Surface* >(sVector,moduleBinUtility);
            // create the layer and push it back
            std::shared_ptr<const Ats::CylinderBounds> cBounds(new Ats::CylinderBounds(layerR, halfZ+layerEnvelopeCoverZ));
            // @TODO overlap descriptor
                                                                                                                                                                  
            // create the layer
            Ats::LayerPtr cLayer = Ats::CylinderLayer::create(nullptr, cBounds, sArray, layerThickness, nullptr, nullptr, Ats::active);
            // push it into the layer vector
            m_cLayers->push_back(cLayer);
        }
    }
    
    // -------------------------------- endcap type layers -----------------------------------------------------------
    // pos/neg layers
    size_t numpnLayers = m_posnegLayerPositionsZ.size();
    if (numpnLayers){
     MSG_DEBUG("Configured to build 2 * " << numpnLayers << " passive positive/negative side layers.");
     m_pLayers = new Ats::LayerVector;
     m_pLayers->reserve(numpnLayers);
     m_nLayers = new Ats::LayerVector;
     m_nLayers->reserve(numpnLayers);
     // loop through
     for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl){
         // some screen output
         MSG_VERBOSE("- build layers " << (2*ipnl) << " and "<<  (2*ipnl)+1 << " at +/- z = " << m_posnegLayerPositionsZ[ipnl]);
         // layer position update
         double layerPosZ                          = m_posnegLayerPositionsZ[ipnl];
         double layerEnvelopeR                     = m_posnegLayerEnvelopeR[ipnl];
         double layerStaggerR                      = m_posnegModuleStaggerR[ipnl];
         // get the min/maximum R
         double layerRmin                          = 10e10;
         double layerRmax                          = 0.;
         // module positioning update
         std::vector<double> layerModuleRadii      = m_posnegModulesRadii[ipnl];
         std::vector<double> layerModulePhi        = m_posnegModulesPhi[ipnl];
         std::vector<double> layerModulePhiStagger = m_posnegMoudleStaggerPhi[ipnl];
         // module description
         std::vector<double> layerModuleMinHalfX   = m_posnegModuleMinHalfX[ipnl];
         std::vector<double> layerModuleMaxHalfX   = m_posnegModuleMaxHalfX[ipnl];
         std::vector<double> layerModuleHalfY      = m_posnegModuleHalfY[ipnl];
         std::vector<double> layerModuleThickness  = m_posnegModuleThickness[ipnl];

         // prepare for the r binning
         std::vector< Ats::SurfaceArray* > pRadialSurfaceArrays;
         std::vector< Ats::SurfaceArray* > nRadialSurfaceArrays;
         std::vector<double> radialBoundariesLow;
         std::vector<double> radialBoudnariesHigh;

         // loop over bins in R
         size_t imodule = 0;
         // layer thickness
         double layerZmax = 0.;
         double layerZmin = 10e10;
         // staggering sterring
         bool rstagger = true;
         // loop over rings
         for (size_t ipnR = 0; ipnR < layerModuleRadii.size(); ++ipnR){
             // incremement @TODO create valid identifier using identifier service
             ++imodule;
             // the actual layer radius & properties of this ring
             double moduleR = layerModuleRadii[ipnR];
             // figure out the staggering
             double moduleZ = layerPosZ;
             moduleZ += rstagger ? 0.5*layerStaggerR : -0.5*layerStaggerR; rstagger = !rstagger;
             // and the bounds
             double moduleMinHalfX  = layerModuleMinHalfX.size() ? layerModuleMinHalfX[ipnR] : 0.;
             double moduleMaxHalfX  = layerModuleMaxHalfX[ipnR];
             double moduleHalfY     = layerModuleHalfY[ipnR];
             double moduleThickness = layerModuleThickness[ipnR];
             // create the bounds
             Ats::PlanarBounds* pBounds =  nullptr;
             if (layerModuleMinHalfX.size())
                   pBounds = new Ats::TrapezoidBounds(moduleMinHalfX, moduleMaxHalfX, moduleHalfY);
             else 
                   pBounds = new Ats::RectangleBounds(moduleMaxHalfX, moduleHalfY); 
             // now create the shared bounds from it
             std::shared_ptr<const Ats::PlanarBounds> moduleBounds(pBounds);
             // stagger in phi
             bool phistagger = true;
             double modulePhiStagger = layerModulePhiStagger[ipnR];
             // phiMin / phiMax
             double minPhi = 10.;
             double maxPhi = -10.;
             // create the
             std::vector< SurfacePosition> nsVector;
             std::vector< SurfacePosition> psVector;
             // now loo over phi
             for (auto& modulePhi : layerModulePhi){
                 // bigger smaller trick on phi
                 takeSmallerBigger(minPhi,maxPhi,modulePhi);
                 // update the module z position
                 moduleZ += phistagger ? 0.5*modulePhiStagger : -0.5*modulePhiStagger; phistagger = !phistagger;
                 // for the z binning
                 takeSmaller(layerZmin, moduleZ-moduleThickness);
                 takeBigger(layerZmax, moduleZ+moduleThickness);
                 // for the disc bounds
                 takeSmaller(layerRmin, moduleR-moduleHalfY);
                 takeBigger(layerRmax, moduleR+moduleHalfY);
                 // the center position of the modules
                 Ats::Vector3D pModuleCenter(moduleR*cos(modulePhi),moduleR*sin(modulePhi),moduleZ);
                 Ats::Vector3D nModuleCenter(moduleR*cos(modulePhi),moduleR*sin(modulePhi),-moduleZ);
                 // the rotation matrix of the module
                 Ats::Vector3D moduleLocalY(cos(modulePhi),sin(modulePhi),0.);
                 Ats::Vector3D pModuleLocalZ(0.,0.,1.); // take different axis to have the same readout direction
                 Ats::Vector3D nModuleLocalZ(0.,0.,-1.); // take different axis to have the same readout direction
                 Ats::Vector3D nModuleLocalX = moduleLocalY.cross(nModuleLocalZ);
                 Ats::Vector3D pModuleLocalX = moduleLocalY.cross(pModuleLocalZ);
                 // local rotation matrices
                 // create the RotationMatrix - negative side
                 Ats::RotationMatrix3D nModuleRotation;
                 nModuleRotation.col(0) = nModuleLocalX;
                 nModuleRotation.col(1) = moduleLocalY;
                 nModuleRotation.col(2) = nModuleLocalZ;
                 // create the RotationMatrix - positive side
                 Ats::RotationMatrix3D pModuleRotation;
                 pModuleRotation.col(0) = pModuleLocalX;
                 pModuleRotation.col(1) = moduleLocalY;
                 pModuleRotation.col(2) = pModuleLocalZ;
                 // the transforms for the two modules
                 std::shared_ptr<Ats::Transform3D> nModuleTransform(new Ats::Transform3D(Ats::getTransformFromRotTransl(nModuleRotation,nModuleCenter)));
                 std::shared_ptr<Ats::Transform3D> pModuleTransform(new Ats::Transform3D(Ats::getTransformFromRotTransl(pModuleRotation,pModuleCenter)));
                 // create the modules identifier @TODO Idenfier service 
                 Identifier nModuleIdentifier(2*imodule);
                 Identifier pModuleIdentifier(2*imodule+1);
                 // create the module 
                 DetectorElement* nmodule = new DetectorElement(nModuleIdentifier, nModuleTransform, moduleBounds, moduleThickness);
                 DetectorElement* pmodule = new DetectorElement(pModuleIdentifier, pModuleTransform, moduleBounds, moduleThickness);
                 
                 // memory management - we need a detector store to hold them somewhere @TODO add detector store facility
                 m_posnegModules.push_back(nmodule);
                 m_posnegModules.push_back(pmodule);
                 
                 // create the surface 
                 nsVector.push_back(SurfacePosition(&nmodule->surface(), nModuleCenter));
                 psVector.push_back(SurfacePosition(&pmodule->surface(), pModuleCenter));                    
             } 
             // create the phi binned array
             double phiStep = (maxPhi-minPhi)/(layerModulePhi.size()-1);
             minPhi -= 0.5*phiStep;
             maxPhi += 0.5*phiStep;
             
             // BinUtilities
             Ats::BinUtility* nphiBinUtility = new Ats::BinUtility(layerModulePhi.size(),minPhi,maxPhi,Ats::closed,Ats::binPhi);
             Ats::BinUtility* pphiBinUtility = new Ats::BinUtility(layerModulePhi.size(),minPhi,maxPhi,Ats::closed,Ats::binPhi);
             nRadialSurfaceArrays.push_back( new Ats::BinnedArray1D<const Ats::Surface* >(nsVector, nphiBinUtility)  );
             pRadialSurfaceArrays.push_back( new Ats::BinnedArray1D<const Ats::Surface* >(psVector, pphiBinUtility)  );
             
         }
         
         // @TODO OverlapDescriptor
         // create the SurfaceArrays
         Ats::SurfaceArray* nSurfaceArray = nullptr;
         Ats::SurfaceArray* pSurfaceArray = nullptr;
         if (nRadialSurfaceArrays.size() == 1 && pRadialSurfaceArrays.size() == 1){
             // just take the one you have
             nSurfaceArray = nRadialSurfaceArrays[0];
             pSurfaceArray = pRadialSurfaceArrays[0];             
         } else {
             // r boundaries ----------------------------------------------
             std::vector< float > rBoundaries = { (float)layerRmin, (float)layerRmax };
             // 
             std::vector< std::pair< Ats::SurfaceArray*, Ats::Vector3D > > pSurfaceArraysPosition;
             std::vector< std::pair< Ats::SurfaceArray*, Ats::Vector3D > > nSurfaceArraysPosition;
             // loop to adjust boundaries
             double innerR      = 0.;
             double outerR      = 0.;
             for (size_t irb = 0; irb < layerModuleRadii.size()-1; ++irb){
                 // needed for boundary calculations
                 innerR      = layerModuleRadii[irb];
                 outerR      = layerModuleRadii[irb+1];
                 double innerHalfY  = layerModuleHalfY[irb];
                 double outerHalfY  = layerModuleHalfY[irb+1];
                 double boundaryR   = 0.5*(innerR+innerHalfY+outerR-outerHalfY);
                 rBoundaries.push_back(boundaryR);
                 // 
                 pSurfaceArraysPosition.push_back( std::pair< Ats::SurfaceArray*, Ats::Vector3D>(pRadialSurfaceArrays[irb], Ats::Vector3D(0.,0.,innerR)) );
                 nSurfaceArraysPosition.push_back( std::pair< Ats::SurfaceArray*, Ats::Vector3D>(nRadialSurfaceArrays[irb], Ats::Vector3D(0.,0.,innerR)) );
             }
             // and the last one
             pSurfaceArraysPosition.push_back( std::pair< Ats::SurfaceArray*, Ats::Vector3D>(pRadialSurfaceArrays[layerModuleRadii.size()-1], Ats::Vector3D(0.,0.,outerR)) );
             nSurfaceArraysPosition.push_back( std::pair< Ats::SurfaceArray*, Ats::Vector3D>(nRadialSurfaceArrays[layerModuleRadii.size()-1], Ats::Vector3D(0.,0.,outerR)) );
             // sort the rBoundaries for the steering bin 
             std::sort(rBoundaries.begin(), rBoundaries.end());
             // building the 1D 1D 
             Ats::BinUtility* nrBinUtility = new Ats::BinUtility(rBoundaries, Ats::open, Ats::binR);
             Ats::BinUtility* prBinUtility = new Ats::BinUtility(rBoundaries, Ats::open, Ats::binR);
             // now create the surface arrays
             nSurfaceArray = new Ats::BinnedArrayArray< const Ats::Surface* >(nSurfaceArraysPosition ,nrBinUtility);
             pSurfaceArray = new Ats::BinnedArrayArray< const Ats::Surface* >(pSurfaceArraysPosition ,prBinUtility);                         
         }
         
         // create the share disc bounds
         std::shared_ptr<const Ats::DiscBounds> dBounds(new Ats::RadialBounds(layerRmin-layerEnvelopeR,layerRmax+layerEnvelopeR));
         
         // layer thickness
         double layerThickness = layerZmax-layerZmin;
         // create the layer transforms
         Ats::Transform3D* nLayerTransform = new Ats::Transform3D(Ats::Transform3D::Identity());
         nLayerTransform->translation() = Ats::Vector3D(0.,0.,-layerPosZ);
         Ats::Transform3D* pLayerTransform = new Ats::Transform3D(Ats::Transform3D::Identity());
         pLayerTransform->translation() = Ats::Vector3D(0.,0.,layerPosZ);

         // create the layers
         Ats::LayerPtr nLayer = Ats::DiscLayer::create(std::shared_ptr<Ats::Transform3D>(nLayerTransform), 
                                                       dBounds,
                                                       nSurfaceArray,
                                                       layerThickness,
                                                       nullptr,
                                                       nullptr,
                                                       Ats::active);
         Ats::LayerPtr pLayer = Ats::DiscLayer::create(std::shared_ptr<Ats::Transform3D>(pLayerTransform), 
                                                       dBounds,
                                                       pSurfaceArray,
                                                       layerThickness,
                                                       nullptr,
                                                       nullptr,
                                                       Ats::active);
         // push it into the layer vector
         m_nLayers->push_back(nLayer);
         m_pLayers->push_back(pLayer);
     }
    }
    
    // everything was successful - let's return back
    return StatusCode::SUCCESS;
}
