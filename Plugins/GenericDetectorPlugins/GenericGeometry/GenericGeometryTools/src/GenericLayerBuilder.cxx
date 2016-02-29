///////////////////////////////////////////////////////////////////
// GenericLayerBuilder.cxx, ATS project
///////////////////////////////////////////////////////////////////

// Core module
#include "Algebra/AlgebraDefinitions.h"
#include "Algebra/AlgebraHelper.h"
// Geometry module
#include "GeometryUtils/BinUtility.h"
#include "GeometryUtils/BinnedArray2D.h"
#include "GeometryUtils/BinnedArray1D.h"
#include "GeometryUtils/BinnedArrayArray.h"
#include "GeometryUtils/ApproachDescriptor.h"
#include "Detector/CylinderLayer.h"
#include "Detector/DiscLayer.h"
#include "Surfaces/RadialBounds.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/PlanarBounds.h"
#include "Surfaces/RectangleBounds.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Material/Material.h"
#include "Material/MaterialProperties.h"
#include "Material/HomogeneousSurfaceMaterial.h"
// A generic detector
#include "GenericGeometryTools/GenericLayerBuilder.h"
#include "GenericDetectorElement/DetectorElement.h"

DECLARE_COMPONENT(Ats::GenericLayerBuilder)

// constructor
Ats::GenericLayerBuilder::GenericLayerBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  Ats::AlgToolBase(t,n,p),
  m_layerIdentification(n),
  m_nLayers(nullptr),
  m_cLayers(nullptr),    
  m_pLayers(nullptr)   
{
    declareInterface<ILayerBuilder>(this);
   
    //  layer identificaiton
    declareProperty("LayerIdentification",                m_layerIdentification);
   
    // the central layers 
    declareProperty("CentralLayerRadii",                  m_centralLayerRadii);
    declareProperty("CentralLayerEnvelopeR",              m_centralLayerEnvelopeR);
    declareProperty("CentralLayerEnvelopeZ",              m_centralLayerEnvelopeZ);
    declareProperty("CentralLayerMaterialConcentration",  m_centralLayerMaterialConcentration);
    declareProperty("CentralLayerMaterialProperties",     m_centralLayerMaterialProperties);    
    declareProperty("CentralLayerModulesPositionPhi",     m_centralModulesPositionPhi);
    declareProperty("CentralLayerMoudlesTiltPhi",         m_centralModulesTiltPhi);        
    declareProperty("CentralLayerModulesPositionZ",       m_centralModulesPositionZ);
    declareProperty("CentralLayerModuleStaggerZ",         m_centralModulesStaggerZ);
    declareProperty("CentralLayerModulesHalfX",           m_centralModuleHalfX);
    declareProperty("CentralLayerModulesHalfY",           m_centralModuleHalfY);
    declareProperty("CentralLayerModulesThickness",       m_centralModuleThickness);
    declareProperty("CentralLayerModulesMaterial",        m_centralModuleMaterial);    
    declareProperty("CentralLayerModulesFrontsideStereo", m_centralModuleFrontsideStereo);
    declareProperty("CentralLayerModulesBacksideStereo",  m_centralModuleBacksideStereo);
    declareProperty("CentralLayerModulesBacksideGap",     m_centralModuleBacksideGap);
    
    // the layers at p/e side 
    declareProperty("PosNegLayerPositionZ",               m_posnegLayerPositionsZ);
    declareProperty("PosNegLayerEnvelopeR",               m_posnegLayerEnvelopeR);
    declareProperty("PosNeglLayerMaterialConcentration",  m_posnegLayerMaterialConcentration);
    declareProperty("PosNeglLayerMaterialProperties",     m_posnegLayerMaterialProperties);    
    declareProperty("PosNegLayerModulesRadii",            m_posnegModulesRadii);
    declareProperty("PosNegLayerModuleStaggerR",          m_posnegModuleStaggerR); 
    declareProperty("PosNegLayerModulesInPhi",            m_posnegModulesInPhi);       
    declareProperty("PosNegLayerModulesPositionPhi",      m_posnegModulesPositionPhiStream);
    declareProperty("PosNegLayerModulesStaggerPhi",       m_posnegMoudleStaggerPhi);
    declareProperty("PosNegLayerModulesMinHalfX",         m_posnegModuleMinHalfX);
    declareProperty("PosNegLayerModulesMaxHalfX",         m_posnegModuleMaxHalfX);
    declareProperty("PosNegLayerModulesHalfY",            m_posnegModuleHalfY);
    declareProperty("PosNegLayerModulesThickness",        m_posnegModuleThickness);
    
}

// destructor
Ats::GenericLayerBuilder::~GenericLayerBuilder()
{}

// initialize
StatusCode Ats::GenericLayerBuilder::initialize()
{
    MSG_DEBUG( "initialize()" );
    
    m_posnegModulesPositionPhi.reserve(m_posnegLayerPositionsZ.size());
    for (size_t idisc = 0; idisc < m_posnegLayerPositionsZ.size(); ++idisc){
         std::vector< std::vector< double > > discPhiPositions;
         discPhiPositions.reserve(m_posnegModulesRadii[idisc].size());
         size_t iphistream = 0;
         for (size_t iring = 0; iring < m_posnegModulesRadii[idisc].size(); ++iring){
             std::vector< double > ringPhiPositions;
             size_t nphiModules = size_t(m_posnegModulesInPhi[idisc][iring]);
             for (size_t iphi = 0; iphi < nphiModules; ++iphi, ++iphistream){
                 ringPhiPositions.push_back(m_posnegModulesPositionPhiStream[idisc][iphistream]);
             }
             discPhiPositions.push_back(ringPhiPositions);    
         }
         m_posnegModulesPositionPhi.push_back(discPhiPositions);
     }
     return constructLayers();
}

//finalize
StatusCode Ats::GenericLayerBuilder::finalize()
{
    MSG_DEBUG( "finalize()" );
    return StatusCode::SUCCESS;
}


StatusCode Ats::GenericLayerBuilder::constructLayers() 
{
    
    // -------------------------------- central layers -----------------------------------------------------------
    typedef std::pair<const Ats::Surface*, Ats::Vector3D> SurfacePosition;
    // the central layers
   size_t numcLayers = m_centralLayerRadii.size();
   if (numcLayers){
       MSG_DEBUG("Configured to build " << numcLayers << " active central layers.");
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
           std::vector<double> phiValues = m_centralModulesPositionPhi[icl];
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
           // temporary cache for neighbor setting
           std::vector< std::vector < DetectorElement* > > neighbourCache;
           std::vector< std::vector < DetectorElement* > > neighbourCacheBackside;    
           // create the Module material from input
           std::shared_ptr<const Ats::SurfaceMaterial> moduleMaterialPtr = nullptr;
           if (m_centralModuleMaterial.size()){
               // get the sensor material - it has to be vectors of 5
               double x0  = m_centralModuleMaterial[icl][0];
               double l0  = m_centralModuleMaterial[icl][1];
               double a   = m_centralModuleMaterial[icl][2];
               double z   = m_centralModuleMaterial[icl][3];
               double rho = m_centralModuleMaterial[icl][4];
               // the moduel moaterial from input 
               Ats::Material moduleMaterial(x0,l0,a,z,rho);
               Ats::MaterialProperties moduleMaterialProperties(moduleMaterial,moduleThickness);
               moduleMaterialPtr = std::shared_ptr<const Ats::SurfaceMaterial>(new Ats::HomogeneousSurfaceMaterial(moduleMaterialProperties));   
           }
           // loop over z module and phi module position
           for (auto& moduleZ : zValues){
               // create the phi neighbours
               std::vector< DetectorElement*> neighbourCachePhi; 
               std::vector< DetectorElement*> neighbourCachePhiBackside; 
               neighbourCachePhi.reserve(phiValues.size());
               neighbourCachePhiBackside.reserve(phiValues.size());
               // create the half length in z
               halfZ = halfZ > moduleZ+moduleHalfY ? halfZ :  moduleZ+moduleHalfY;
               // loop of phi values
               for (auto& modulePhi : phiValues){
                   // min/max phi
                   takeSmallerBigger(minPhi, maxPhi, modulePhi);

                   // stagger the modules
                   double moduleR = layerR;
                   moduleR += stagger ? -0.5*layerModleStaggerZ : 0.5 * layerModleStaggerZ; stagger = !stagger;
                   // the position of the module
                   Ats::Vector3D moduleCenter(moduleR*cos(modulePhi), moduleR*sin(modulePhi), moduleZ);  
                   // normal vectorof the surface
                   Ats::Vector3D moduleLocalZ(cos(modulePhi+modulePhiTilt),sin(modulePhi+modulePhiTilt), 0.);
                   Ats::Vector3D moduleLocalY(0.,0.,1);
                   Ats::Vector3D moduleLocalX(-sin(modulePhi+modulePhiTilt),cos(modulePhi+modulePhiTilt),0.);
                   // create the RotationMatrix
                   Ats::RotationMatrix3D moduleRotation;
                   moduleRotation.col(0) = moduleLocalX;
                   moduleRotation.col(1) = moduleLocalY;
                   moduleRotation.col(2) = moduleLocalZ;
                   // get the moduleTransform
                   std::shared_ptr<Ats::Transform3D> moduleTransform(new Ats::Transform3D(Ats::getTransformFromRotTransl(moduleRotation,moduleCenter)));
                   // stereo angle applied
                   if (m_centralModuleFrontsideStereo.size() && m_centralModuleFrontsideStereo[icl] != 0.){
                       // twist by the stereo angle
                       double stereo = m_centralModuleFrontsideStereo[icl];
                       (*moduleTransform.get()) *= Ats::AngleAxis3D(-stereo, Ats::Vector3D::UnitZ());
                   }
                   // count the modules
                   ++imodule;
                   Identifier moduleIdentifier = Identifier(Identifier::value_type(imodule));
                   // create the module 
                   DetectorElement* module = new DetectorElement(moduleIdentifier, moduleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                   // radial extend
                   moduleRadialExtend(*module,moduleThickness,moduleHalfX,moduleHalfY,layerMinR,layerMaxR);
                   // register the module to the cache module
                   neighbourCachePhi.push_back(module);
                   // create the surface 
                   sVector.push_back(SurfacePosition(&module->surface(),moduleCenter));
                   // memory management - we need a detector store to hold them somewhere @TODO detector store facility
                   m_centralModules.push_back(module);
                   // and the backside one (if configured to do so)
                   if (m_centralModuleBacksideGap.size()){
                       // ncrease the counter @TODO switch to identifier service
                       ++imodule;
                       // create the module identifier
                       moduleIdentifier = Identifier(Identifier::value_type(imodule));
                       moduleCenter = moduleCenter + m_centralModuleBacksideGap[icl]*moduleLocalZ;                  
                       moduleTransform = std::shared_ptr<Ats::Transform3D>(new Ats::Transform3D(Ats::getTransformFromRotTransl(moduleRotation,moduleCenter)));
                       // apply the stereo
                       if (m_centralModuleBacksideStereo.size()){
                           // twist by the stereo angle
                           double stereoBackSide = m_centralModuleBacksideStereo[icl];
                           (*moduleTransform.get()) *= Ats::AngleAxis3D(-stereoBackSide, Ats::Vector3D::UnitZ());
                       }
                       // everything is set for the next module
                       DetectorElement* bsmodule = new DetectorElement(moduleIdentifier, moduleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                       // radial extend
                       moduleRadialExtend(*bsmodule,moduleThickness,moduleHalfX,moduleHalfY,layerMinR,layerMaxR);
                       // register the module to the cache module
                       neighbourCachePhiBackside.push_back(bsmodule);
                       // memory management - we need a detector store to hold them somewhere @TODO detector store facility
                       m_centralModules.push_back(bsmodule);
                   }
               }
               // register the phi vector to the module cache
               neighbourCache.push_back(neighbourCachePhi);
               if (neighbourCachePhiBackside.size()) 
                   neighbourCacheBackside.push_back(neighbourCachePhiBackside);

           }
           // create the neighbor Cache
           MSG_VERBOSE("Filling the neighbour cache for the overlap descriptor.");
           // --- interlude --- fill the neighbour cache
           // number of z modules
           size_t zModules = neighbourCache.size();
           for (size_t iz = 0; iz < zModules; ++iz){
               // number of phi modules
               size_t phiModules = neighbourCache[iz].size();
               // the previous z and next z
               int pz = (int(iz) > 0) ? iz-1 : -1;
               int nz = (int(iz) < int(zModules-1)) ? iz+1 : -1; 
               // now loop over the phi values
               for (size_t iphi = 0; iphi < phiModules; ++iphi){
                   // create the neighbours
                   std::vector<const Ats::DetectorElementBase*> neighbours;
                   // closed for phi
                   size_t nphi = (iphi < (phiModules-1)) ? iphi+1 : 0;
                   size_t pphi = (iphi > 0) ? iphi-1 : phiModules-1;
                   // fill neighbours
                   // get the current z with next in phi
                   neighbours.push_back(neighbourCache[iz][pphi]);
                   neighbours.push_back(neighbourCache[iz][nphi]);
                   // get the z neighbours
                   if (pz > 0){
                       neighbours.push_back(neighbourCache[pz][pphi]);
                       neighbours.push_back(neighbourCache[pz][iphi]);
                       neighbours.push_back(neighbourCache[pz][nphi]);
                   }
                   if (nz > 0){
                       neighbours.push_back(neighbourCache[nz][pphi]);
                       neighbours.push_back(neighbourCache[nz][iphi]);
                       neighbours.push_back(neighbourCache[nz][nphi]);
                   }                    
                   // fill backside neighbours if created
                   if (neighbourCacheBackside.size()){
                       neighbours.push_back(neighbourCacheBackside[iz][pphi]);
                       neighbours.push_back(neighbourCacheBackside[iz][iphi]); 
                       neighbours.push_back(neighbourCacheBackside[iz][nphi]);
                       // get the z neighbours
                       if (pz > 0){
                           neighbours.push_back(neighbourCacheBackside[pz][pphi]);
                           neighbours.push_back(neighbourCacheBackside[pz][iphi]);
                           neighbours.push_back(neighbourCacheBackside[pz][nphi]);
                       }
                       if (nz > 0){
                           neighbours.push_back(neighbourCacheBackside[nz][pphi]);
                           neighbours.push_back(neighbourCacheBackside[nz][iphi]);
                           neighbours.push_back(neighbourCacheBackside[nz][nphi]);
                       } 
                   }
                   // now register it to the volume
                   if (neighbourCache[iz][iphi]) neighbourCache[iz][iphi]->registerNeighbours(neighbours);                    
               }
           }

           // harmonize the phi boundaries 
           double phiStep = (maxPhi-minPhi)/(phiValues.size()-1);
           minPhi -= 0.5*phiStep;
           maxPhi += 0.5*phiStep;
           // layer thickness
           double layerThickness = (layerMaxR-layerMinR);
           MSG_VERBOSE("-        with z min/max = " << -halfZ << " / " << halfZ);
           MSG_VERBOSE("-        with R min/max = " << layerMinR << " / " << layerMaxR);
           MSG_VERBOSE("- and number of modules = " << sVector.size() << " ( " << phiValues.size() << " x " << zValues.size() << ")");

           // create the binUtility
           Ats::BinUtility* moduleBinUtility = new Ats::BinUtility(phiValues.size(), minPhi, maxPhi, Ats::closed, Ats::binPhi);
           (*moduleBinUtility) += Ats::BinUtility(zValues.size(), -halfZ, halfZ, Ats::open, Ats::binZ);
           // create the surface array 
           sArray = new Ats::BinnedArray2D< const Ats::Surface* >(sVector,moduleBinUtility);
           // create the layer and push it back
           std::shared_ptr<const Ats::CylinderBounds> cBounds(new Ats::CylinderBounds(layerR, halfZ+layerEnvelopeCoverZ));
           // create the layer
           Ats::LayerPtr cLayer = Ats::CylinderLayer::create(nullptr, cBounds, sArray, layerThickness, nullptr, nullptr, Ats::active);
           // the layer is built le't see if it needs material
           if (m_centralLayerMaterialProperties.size()){
               // get the material from configuration
               double lMaterialThickness = m_centralLayerMaterialProperties[icl][0];
               double lMaterialX0        = m_centralLayerMaterialProperties[icl][1];
               double lMaterialL0        = m_centralLayerMaterialProperties[icl][2];
               double lMaterialA         = m_centralLayerMaterialProperties[icl][3];
               double lMaterialZ         = m_centralLayerMaterialProperties[icl][4];
               double lMaterialRho       = m_centralLayerMaterialProperties[icl][5];
               Ats::MaterialProperties layerMaterialProperties(lMaterialThickness,lMaterialX0,lMaterialL0,lMaterialA,lMaterialZ,lMaterialRho);
               std::shared_ptr<const Ats::SurfaceMaterial> layerMaterialPtr(new Ats::HomogeneousSurfaceMaterial(layerMaterialProperties));   
               // get the approach descriptor - at this stage we know that the approachDescriptor exists
               auto approachSurfaces = cLayer->approachDescriptor()->containedSurfaces();
               if (m_centralLayerMaterialConcentration[icl] > 0)
                   approachSurfaces[1]->setSurfaceMaterial(layerMaterialPtr);
               else 
                   approachSurfaces[0]->setSurfaceMaterial(layerMaterialPtr);
           }
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
         auto layerModuleRadii                     = m_posnegModulesRadii[ipnl];
         auto layerModulePositionsPhi              = m_posnegModulesPositionPhi[ipnl];
         auto layerModulePhiStagger                = m_posnegMoudleStaggerPhi[ipnl];
         // module description
         auto layerModuleMinHalfX                  = m_posnegModuleMinHalfX[ipnl];
         auto layerModuleMaxHalfX                  = m_posnegModuleMaxHalfX[ipnl];
         auto layerModuleHalfY                     = m_posnegModuleHalfY[ipnl];
         auto layerModuleThickness                 = m_posnegModuleThickness[ipnl];

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
             // the phi module of this ring
             auto ringModulePositionsPhi = layerModulePositionsPhi[ipnR];
             // now loo over phi
             for (auto& modulePhi : ringModulePositionsPhi){
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
                 Identifier nModuleIdentifier = Identifier(Identifier::value_type(2*imodule));
                 Identifier pModuleIdentifier = Identifier(Identifier::value_type(2*imodule+1));
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
             double phiStep = (maxPhi-minPhi)/(ringModulePositionsPhi.size()-1);
             minPhi -= 0.5*phiStep;
             maxPhi += 0.5*phiStep;
             
             // BinUtilities
             Ats::BinUtility* nphiBinUtility = new Ats::BinUtility(ringModulePositionsPhi.size(),minPhi,maxPhi,Ats::closed,Ats::binPhi);
             Ats::BinUtility* pphiBinUtility = new Ats::BinUtility(ringModulePositionsPhi.size(),minPhi,maxPhi,Ats::closed,Ats::binPhi);
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

void Ats::GenericLayerBuilder::moduleRadialExtend(const DetectorElement& detElement, double thickness, double halfX, double halfY, double& rmin, double& rmax){
    
    // brute force method to find the radial extends for the layer construction
    
    Vector3D localX = detElement.transform().rotation().col(0);
    Vector3D localY = detElement.transform().rotation().col(1);
    Vector3D localZ = detElement.transform().rotation().col(2);
    // center
    Vector3D center = detElement.transform().translation();
    
    // construct the edges 
    Vector3D edge000 = center-0.5*thickness*localZ-halfX*localX-halfY*localY;
    Vector3D edge001 = center-0.5*thickness*localZ-halfX*localX+halfY*localY;
    Vector3D edge010 = center-0.5*thickness*localZ+halfX*localX-halfY*localY;
    Vector3D edge011 = center-0.5*thickness*localZ+halfX*localX+halfY*localY;
    Vector3D edge100 = center+0.5*thickness*localZ-halfX*localX-halfY*localY;
    Vector3D edge101 = center+0.5*thickness*localZ-halfX*localX+halfY*localY;
    Vector3D edge110 = center+0.5*thickness*localZ+halfX*localX-halfY*localY;
    Vector3D edge111 = center+0.5*thickness*localZ+halfX*localX+halfY*localY;
    // layer min / max 
    takeSmallerBigger(rmin,rmax,edge000.perp());
    takeSmallerBigger(rmin,rmax,edge001.perp());
    takeSmallerBigger(rmin,rmax,edge010.perp());
    takeSmallerBigger(rmin,rmax,edge011.perp());
    takeSmallerBigger(rmin,rmax,edge100.perp());
    takeSmallerBigger(rmin,rmax,edge101.perp());
    takeSmallerBigger(rmin,rmax,edge110.perp());
    takeSmallerBigger(rmin,rmax,edge111.perp());
}


