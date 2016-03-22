///////////////////////////////////////////////////////////////////
// GenericLayerBuilder.cxx, ACTS project
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
#include "GenericDetectorElement/GenericDetectorElement.h"

DECLARE_TOOL_FACTORY(Acts::GenericLayerBuilder)

// constructor
Acts::GenericLayerBuilder::GenericLayerBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  AlgToolBase(t,n,p),
  m_layerIdentification(n),
  m_nLayers(nullptr),
  m_cLayers(nullptr),    
  m_pLayers(nullptr),
  m_approachSurfaceEnvelope(0.1),
  m_centralPassiveLayerBuilder(""),
  m_posnegPassiveLayerBuilder("")
{
    declareInterface<ILayerBuilder>(this);
   
    //  layer identificaiton
    declareProperty("LayerIdentification",                m_layerIdentification);
   
    // approach surface envelope
    declareProperty("ApproachSurfaceEnvelope",            m_approachSurfaceEnvelope);
   
    // the central layers 
    declareProperty("CentralLayerRadii",                  m_centralLayerRadii);
    declareProperty("CentralLayerEnvelopeR",              m_centralLayerEnvelopeR);
    declareProperty("CentralLayerEnvelopeZ",              m_centralLayerEnvelopeZ);
    declareProperty("CentralLayerMaterialConcentration",  m_centralLayerMaterialConcentration);
    declareProperty("CentralLayerMaterialProperties",     m_centralLayerMaterialProperties);    
    declareProperty("CentralLayerModulesPositionPhi",     m_centralModulePositionPhi);
    declareProperty("CentralLayerMoudlesTiltPhi",         m_centralModuleTiltPhi);        
    declareProperty("CentralLayerModulesPositionZ",       m_centralModulePositionZ);
    declareProperty("CentralLayerModuleStaggerZ",         m_centralModuleStaggerZ);
    declareProperty("CentralLayerModulesHalfX",           m_centralModuleHalfX);
    declareProperty("CentralLayerModulesHalfY",           m_centralModuleHalfY);
    declareProperty("CentralLayerModulesThickness",       m_centralModuleThickness);
    declareProperty("CentralLayerModulesMaterial",        m_centralModuleMaterial);    
    declareProperty("CentralLayerModulesFrontsideStereo", m_centralModuleFrontsideStereo);
    declareProperty("CentralLayerModulesBacksideStereo",  m_centralModuleBacksideStereo);
    declareProperty("CentralLayerModulesBacksideGap",     m_centralModuleBacksideGap);
    declareProperty("CentralPassiveLayerBuilder",         m_centralPassiveLayerBuilder);
    
    // the layers at p/e side 
    declareProperty("PosNegLayerPositionZ",               m_posnegLayerPositionsZ);
    declareProperty("PosNegLayerEnvelopeR",               m_posnegLayerEnvelopeR);
    declareProperty("PosNegLayerMaterialConcentration",   m_posnegLayerMaterialConcentration);
    declareProperty("PosNegLayerMaterialProperties",      m_posnegLayerMaterialProperties);    
    declareProperty("PosNegLayerModulesRadii",            m_posnegModuleRadii);
    declareProperty("PosNegLayerModuleStaggerR",          m_posnegModuleStaggerR); 
    declareProperty("PosNegLayerModulesInPhi",            m_posnegModuleInPhi);       
    declareProperty("PosNegLayerModulesPositionPhi",      m_posnegModulePositionPhiStream);
    declareProperty("PosNegLayerModulesStaggerPhi",       m_posnegMoudleStaggerPhi);
    declareProperty("PosNegLayerModulesMinHalfX",         m_posnegModuleMinHalfX);
    declareProperty("PosNegLayerModulesMaxHalfX",         m_posnegModuleMaxHalfX);
    declareProperty("PosNegLayerModulesHalfY",            m_posnegModuleHalfY);
    declareProperty("PosNegLayerModulesThickness",        m_posnegModuleThickness);
    declareProperty("PosNegLayerModulesMaterial",         m_posnegModuleMaterialStream);
    declareProperty("PosNegModulesFrontsideStereo",       m_posnegModuleFrontsideStereo);
    declareProperty("PosNegModulesBacksideStereo",        m_posnegModuleBacksideStereo);
    declareProperty("PosNegModulesBacksideGap",           m_posnegModuleBacksideGap);
    declareProperty("PosNegPassiveLayerBuilder",          m_posnegPassiveLayerBuilder);
    
}

// destructor
Acts::GenericLayerBuilder::~GenericLayerBuilder()
{}

// initialize
StatusCode Acts::GenericLayerBuilder::initialize()
{
    MSG_DEBUG( "initialize()" );
    //Tool needs to be initialized
    if (!AlgToolBase::initialize()) return StatusCode::FAILURE;
    
    // retrieve the passive layer builders if there
    RETRIEVE_NONEMPTY_FATAL(m_centralPassiveLayerBuilder);
    RETRIEVE_NONEMPTY_FATAL(m_posnegPassiveLayerBuilder);

    // convert the stream into ordered lists - no checking for dimensions
    m_posnegModulePositionPhi.reserve(m_posnegLayerPositionsZ.size());
    m_posnegModuleMaterial.reserve(m_posnegLayerPositionsZ.size());
    for (size_t idisc = 0; idisc < m_posnegLayerPositionsZ.size(); ++idisc){
         // that's for the phi positions 
         std::vector< std::vector< double > > discPhiPositions;      
         discPhiPositions.reserve(m_posnegModuleRadii[idisc].size());
         // that's for the module material
         std::vector < std::vector< double > > discModuleMaterial;
         discModuleMaterial.reserve(m_posnegModuleRadii[idisc].size());
         size_t iphistream = 0;
         size_t imstream   = 0;
         for (size_t iring = 0; iring < m_posnegModuleRadii[idisc].size(); ++iring){
             // resolve the phi positions
             std::vector< double > ringPhiPositions;
             size_t nphiModules = size_t(m_posnegModuleInPhi[idisc][iring]);
             for (size_t iphi = 0; iphi < nphiModules; ++iphi, ++iphistream){
                 ringPhiPositions.push_back(m_posnegModulePositionPhiStream[idisc][iphistream]);
             }
             discPhiPositions.push_back(ringPhiPositions); 
             
             if (m_posnegModuleMaterialStream.size()){
                // resolve the material
                std::vector< double > ringModuleMaterial;
                for (size_t im = 0; im < 5; ++im, ++imstream)
                    ringModuleMaterial.push_back(m_posnegModuleMaterialStream[idisc][imstream]);
                discModuleMaterial.push_back(ringModuleMaterial);
            }   
         }
         m_posnegModulePositionPhi.push_back(discPhiPositions);
         if (discModuleMaterial.size()) 
             m_posnegModuleMaterial.push_back(discModuleMaterial);
     }
     return constructLayers();
}

//finalize
StatusCode Acts::GenericLayerBuilder::finalize()
{
    MSG_DEBUG( "finalize()" );
    return StatusCode::SUCCESS;
}


StatusCode Acts::GenericLayerBuilder::constructLayers() 
{
    
   // -------------------------------- central layers -----------------------------------------------------------
   typedef std::pair<const Surface*, Vector3D> SurfacePosition;
   // the central layers
   size_t numcLayers = m_centralLayerRadii.size();
   if (numcLayers){
       MSG_DEBUG("Configured to build " << numcLayers << " active central layers.");
       m_cLayers = new LayerVector;
       m_cLayers->reserve(numcLayers);
       // loop through
       for (size_t icl = 0; icl < numcLayers; ++icl){
           // layer R/Z
           double layerR      = m_centralLayerRadii[icl];
           double layerHalfZ  = 0.;
           double minPhi      = 10.;
           double maxPhi      = -10.;
           // some screen output
           MSG_VERBOSE("- build layer " << icl << " with target radius = " << layerR);
           // create the modules & surface array 
           SurfaceArray* sArray = nullptr;
           // surface vector 
           std::vector<SurfacePosition> sVector;
           // z/phi values for this layer
           std::vector<double> zValues   = m_centralModulePositionZ[icl];
           std::vector<double> phiValues = m_centralModulePositionPhi[icl];
           std::sort(phiValues.begin(),phiValues.end());
           std::sort(zValues.begin(),zValues.end());
           // envelope stagger & cover
           double layerModleStaggerZ  = m_centralModuleStaggerZ.size() ? m_centralModuleStaggerZ[icl] : 0.;
           double layerEnvelopeCoverZ = m_centralLayerEnvelopeZ.size() ? m_centralLayerEnvelopeZ[icl] : 0.;
           double layerMinR = 10e10;
           double layerMaxR = 0;
           double layerMinZ = 10e10;
           double layerMaxZ = -10e10;
           // module size & tilt
           double modulePhiTilt   = m_centralModuleTiltPhi[icl]; 
           double moduleHalfX     = m_centralModuleHalfX[icl];
           double moduleHalfY     = m_centralModuleHalfY[icl];
           double moduleThickness = m_centralModuleThickness[icl];
           // create the shared module 
           std::shared_ptr<const PlanarBounds> moduleBounds(new RectangleBounds(moduleHalfX,moduleHalfY));
           // now create the modules and surfaces 
           bool stagger = false;
           // Identifier @TODO unique Identifier - use a GenericDetector identifier
           size_t imodule = 0;
           sVector.reserve(zValues.size()*phiValues.size());
           // temporary cache for neighbor setting
           std::vector< std::vector < GenericDetectorElement* > > neighbourCache;
           std::vector< std::vector < GenericDetectorElement* > > neighbourCacheBackside;    
           // create the Module material from input
           std::shared_ptr<const SurfaceMaterial> moduleMaterialPtr = nullptr;
           if (m_centralModuleMaterial.size()){
               // get the sensor material - it has to be vectors of 5
               double x0  = m_centralModuleMaterial[icl][0];
               double l0  = m_centralModuleMaterial[icl][1];
               double a   = m_centralModuleMaterial[icl][2];
               double z   = m_centralModuleMaterial[icl][3];
               double rho = m_centralModuleMaterial[icl][4];
               // the moduel moaterial from input 
               Material moduleMaterial(x0,l0,a,z,rho);
               MaterialProperties moduleMaterialProperties(moduleMaterial,moduleThickness);
               moduleMaterialPtr = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(moduleMaterialProperties));   
           }
           // loop over z module and phi module position
           for (auto& moduleZ : zValues){
               // create the phi neighbours
               std::vector< GenericDetectorElement*> neighbourCachePhi; 
               std::vector< GenericDetectorElement*> neighbourCachePhiBackside; 
               neighbourCachePhi.reserve(phiValues.size());
               neighbourCachePhiBackside.reserve(phiValues.size());
               // loop of phi values
               for (auto& modulePhi : phiValues){
                   // min/max phi
                   takeSmallerBigger(minPhi, maxPhi, modulePhi);
                   // stagger the modules
                   double moduleR = layerR;
                   moduleR += stagger ? -0.5*layerModleStaggerZ : 0.5 * layerModleStaggerZ; stagger = !stagger;
                   // the position of the module
                   Vector3D moduleCenter(moduleR*cos(modulePhi), moduleR*sin(modulePhi), moduleZ);  
                   // normal vectorof the surface
                   Vector3D moduleLocalZ(cos(modulePhi+modulePhiTilt),sin(modulePhi+modulePhiTilt), 0.);
                   Vector3D moduleLocalY(0.,0.,1);
                   Vector3D moduleLocalX(-sin(modulePhi+modulePhiTilt),cos(modulePhi+modulePhiTilt),0.);
                   // create the RotationMatrix
                   RotationMatrix3D moduleRotation;
                   moduleRotation.col(0) = moduleLocalX;
                   moduleRotation.col(1) = moduleLocalY;
                   moduleRotation.col(2) = moduleLocalZ;
                   // get the moduleTransform
                   std::shared_ptr<Transform3D> moduleTransform(new Transform3D(getTransformFromRotTransl(moduleRotation,moduleCenter)));
                   // stereo angle applied
                   if (m_centralModuleFrontsideStereo.size() && m_centralModuleFrontsideStereo[icl] != 0.){
                       // twist by the stereo angle
                       double stereo = m_centralModuleFrontsideStereo[icl];
                       (*moduleTransform.get()) *= AngleAxis3D(-stereo, Vector3D::UnitZ());
                   }
                   // count the modules
                   ++imodule;
                   Identifier moduleIdentifier = Identifier(Identifier::value_type(imodule));
                   // create the module 
                   GenericDetectorElement* module = new GenericDetectorElement(moduleIdentifier, moduleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                   // radial extend
                   moduleExtend(*module,moduleThickness,moduleHalfX,moduleHalfX,moduleHalfY,layerMinR,layerMaxR,layerMinZ,layerMaxZ);
                   // register the module to the cache module
                   neighbourCachePhi.push_back(module);
                   // create the surface 
                   sVector.push_back(SurfacePosition(&module->surface(),moduleCenter));
                   // memory management - we need a detector store to hold them somewhere @TODO detector store facility
                   m_centralModule.push_back(module);
                   // and the backside one (if configured to do so)
                   if (m_centralModuleBacksideGap.size()){
                       // ncrease the counter @TODO switch to identifier service
                       ++imodule;
                       // create the module identifier
                       moduleIdentifier = Identifier(Identifier::value_type(imodule));
                       moduleCenter = moduleCenter + m_centralModuleBacksideGap[icl]*moduleLocalZ;                  
                       moduleTransform = std::shared_ptr<Transform3D>(new Transform3D(getTransformFromRotTransl(moduleRotation,moduleCenter)));
                       // apply the stereo
                       if (m_centralModuleBacksideStereo.size()){
                           // twist by the stereo angle
                           double stereoBackSide = m_centralModuleBacksideStereo[icl];
                           (*moduleTransform.get()) *= AngleAxis3D(-stereoBackSide, Vector3D::UnitZ());
                       }
                       // everything is set for the next module
                       GenericDetectorElement* bsmodule = new GenericDetectorElement(moduleIdentifier, moduleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                       // radial extend
                       moduleExtend(*bsmodule,moduleThickness,moduleHalfX,moduleHalfX,moduleHalfY,layerMinR,layerMaxR,layerMinZ,layerMaxZ);
                       // register the module to the cache module
                       neighbourCachePhiBackside.push_back(bsmodule);
                       // memory management - we need a detector store to hold them somewhere @TODO detector store facility
                       m_centralModule.push_back(bsmodule);
                   }
               }
               // register the phi vector to the module cache
               neighbourCache.push_back(neighbourCachePhi);
               if (neighbourCachePhiBackside.size()) 
                   neighbourCacheBackside.push_back(neighbourCachePhiBackside);

           }
           // create the neighbor Cache
           MSG_VERBOSE("- filling the neighbour cache for the overlap descriptor.");
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
                   std::vector<const DetectorElementBase*> neighbours;
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
                   if (neighbourCache[iz][iphi]) {
                       neighbourCache[iz][iphi]->registerNeighbours(neighbours);
                       // now register the almost same stuff to the backside module
                       if (neighbourCacheBackside.size() && neighbourCacheBackside[iz][iphi]){
                           // these are the neighbours of the backside module
                           std::vector<const DetectorElementBase*> bsneighbours; 
                           bsneighbours.reserve(neighbours.size());
                           for (auto& nmodule : neighbours){
                               if (nmodule != neighbourCacheBackside[iz][iphi])
                                   bsneighbours.push_back(nmodule);
                           }
                           // now register the final one
                           bsneighbours.push_back(neighbourCache[iz][iphi]);
                           neighbourCacheBackside[iz][iphi]->registerNeighbours(bsneighbours);
                       }
                   }
               }
           }

           // harmonize the phi boundaries 
           double phiStep = (maxPhi-minPhi)/(phiValues.size()-1);
           minPhi -= 0.5*phiStep;
           maxPhi += 0.5*phiStep;
           // layer thickness
           double layerThickness = (layerMaxR-layerMinR)+2*m_approachSurfaceEnvelope;
           // layer half z 
           takeBigger(layerHalfZ,layerMaxZ); 
           // adjust the layer radius 
           layerR = 0.5*(layerMinR+layerMaxR);
           MSG_VERBOSE("-        with layer R   = " << layerR);
           MSG_VERBOSE("-        from R min/max = " << layerMinR << " / " << layerMaxR);
           MSG_VERBOSE("-        with z min/max = " << -layerHalfZ << " / " << layerHalfZ);
           MSG_VERBOSE("-        # of modules   = " << sVector.size() << " ( " << phiValues.size() << " x " << zValues.size() << ")");

           // create the binUtility
           BinUtility* moduleBinUtility = new BinUtility(phiValues.size(), minPhi, maxPhi, closed, binPhi);
           (*moduleBinUtility) += BinUtility(zValues.size(), -layerHalfZ, layerHalfZ, open, binZ);
           // create the surface array 
           sArray = new BinnedArray2D< const Surface* >(sVector,moduleBinUtility);
           // create the layer and push it back
           std::shared_ptr<const CylinderBounds> cBounds(new CylinderBounds(layerR, layerHalfZ+layerEnvelopeCoverZ));
           // create the layer
           LayerPtr cLayer = CylinderLayer::create(nullptr, cBounds, sArray, layerThickness, nullptr, nullptr, active);
           // the layer is built le't see if it needs material
           if (m_centralLayerMaterialProperties.size()){
               // get the material from configuration
               double lMaterialThickness = m_centralLayerMaterialProperties[icl][0];
               double lMaterialX0        = m_centralLayerMaterialProperties[icl][1];
               double lMaterialL0        = m_centralLayerMaterialProperties[icl][2];
               double lMaterialA         = m_centralLayerMaterialProperties[icl][3];
               double lMaterialZ         = m_centralLayerMaterialProperties[icl][4];
               double lMaterialRho       = m_centralLayerMaterialProperties[icl][5];
               MaterialProperties layerMaterialProperties(lMaterialThickness,lMaterialX0,lMaterialL0,lMaterialA,lMaterialZ,lMaterialRho);
               std::shared_ptr<const SurfaceMaterial> layerMaterialPtr(new HomogeneousSurfaceMaterial(layerMaterialProperties));   
               // central material
               if (m_centralLayerMaterialConcentration[icl] == 0.){
                   // the layer surface is the material surface
                   MSG_VERBOSE("- and material at central surface at radius =  " << cLayer->surfaceRepresentation().bounds().r() );
                   cLayer->surfaceRepresentation().setSurfaceMaterial(layerMaterialPtr);
               } else {
                   // approach surface material
                   // get the approach descriptor - at this stage we know that the approachDescriptor exists
                   auto approachSurfaces = cLayer->approachDescriptor()->containedSurfaces();
                   if (m_centralLayerMaterialConcentration[icl] > 0){
                       approachSurfaces[1]->setSurfaceMaterial(layerMaterialPtr);
                       MSG_VERBOSE("- and material at outer approach surfaceat radius =  " << approachSurfaces[1]->bounds().r() );
                   }
                   else {
                       approachSurfaces[0]->setSurfaceMaterial(layerMaterialPtr);
                       MSG_VERBOSE("- and material at inner approach surface at radius =  " << approachSurfaces[0]->bounds().r() );
                   }
               }
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
     m_pLayers = new LayerVector;
     m_pLayers->reserve(numpnLayers);
     m_nLayers = new LayerVector;
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
         auto layerModuleRadii                     = m_posnegModuleRadii[ipnl];
         auto layerModulePositionsPhi              = m_posnegModulePositionPhi[ipnl];
         auto layerModulePhiStagger                = m_posnegMoudleStaggerPhi[ipnl];
         // module description
         auto layerModuleMinHalfX                  = m_posnegModuleMinHalfX[ipnl];
         auto layerModuleMaxHalfX                  = m_posnegModuleMaxHalfX[ipnl];
         auto layerModuleHalfY                     = m_posnegModuleHalfY[ipnl];
         auto layerModuleThickness                 = m_posnegModuleThickness[ipnl];

         // prepare for the r binning
         std::vector< SurfaceArray* > pRadialSurfaceArrays;
         std::vector< SurfaceArray* > nRadialSurfaceArrays;
         std::vector<double> radialBoundariesLow;
         std::vector<double> radialBoudnariesHigh;

         // the phiR Cache              
         std::vector< std::vector< GenericDetectorElement*> > nneighbourCache; 
         std::vector< std::vector< GenericDetectorElement*> > nneighbourCacheBackside; 
         std::vector< std::vector< GenericDetectorElement*> > pneighbourCache; 
         std::vector< std::vector< GenericDetectorElement*> > pneighbourCacheBackside; 

         // loop over bins in R
         size_t imodule = 0;
         // layer thickness
         double layerZmax = 0.;
         double layerZmin = 10e10;
         // staggering sterring
         bool rstagger = true;
         // screen output
         MSG_VERBOSE("This pair of discs has " << layerModuleRadii.size() << " rings.");
         // loop over rings
         for (size_t ipnR = 0; ipnR < layerModuleRadii.size(); ++ipnR){

             // the phi Cache              
             std::vector< GenericDetectorElement*> nneighbourCachePhi; 
             std::vector< GenericDetectorElement*> nneighbourCachePhiBackside; 
             std::vector< GenericDetectorElement*> pneighbourCachePhi; 
             std::vector< GenericDetectorElement*> pneighbourCachePhiBackside; 

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
             
             MSG_VERBOSE("Ring - " << ipnR << " - Checking for sensor material to be built for sensors");
             // create the Module material from input
             std::shared_ptr<const SurfaceMaterial> moduleMaterialPtr = nullptr;
             if (m_posnegModuleMaterial.size()){
                 // get the sensor material - it has to be vectors of 5
                 double x0  = m_posnegModuleMaterial[ipnl][ipnR][0];
                 double l0  = m_posnegModuleMaterial[ipnl][ipnR][1];
                 double a   = m_posnegModuleMaterial[ipnl][ipnR][2];
                 double z   = m_posnegModuleMaterial[ipnl][ipnR][3];
                 double rho = m_posnegModuleMaterial[ipnl][ipnR][4];
                 // the moduel moaterial from input 
                 Material moduleMaterial(x0,l0,a,z,rho);
                 MaterialProperties moduleMaterialProperties(moduleMaterial,moduleThickness);
                 // and create the shared pointer
                 moduleMaterialPtr = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(moduleMaterialProperties));   
             }
             
             // create the bounds
             PlanarBounds* pBounds =  nullptr;
             if (layerModuleMinHalfX.size())
                   pBounds = new TrapezoidBounds(moduleMinHalfX, moduleMaxHalfX, moduleHalfY);
             else 
                   pBounds = new RectangleBounds(moduleMaxHalfX, moduleHalfY); 
             // now create the shared bounds from it
             std::shared_ptr<const PlanarBounds> moduleBounds(pBounds);
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
             
             MSG_VERBOSE("Ring - " << ipnR << " - has " << ringModulePositionsPhi.size() << " phi modules.");
             
             // now loo over phi
             for (auto& modulePhi : ringModulePositionsPhi){
                 // bigger smaller trick on phi
                 takeSmallerBigger(minPhi,maxPhi,modulePhi);
                 // update the module z position
                 moduleZ += phistagger ? 0.5*modulePhiStagger : -0.5*modulePhiStagger; phistagger = !phistagger;
                 // for the z binning
                 takeSmaller(layerZmin, moduleZ-0.5*moduleThickness);
                 takeBigger(layerZmax, moduleZ+0.5*moduleThickness);
                 // for the disc bounds
                 takeSmaller(layerRmin, moduleR-moduleHalfY);
                 takeBigger(layerRmax, moduleR+moduleHalfY);
                 // the center position of the modules
                 Vector3D pModuleCenter(moduleR*cos(modulePhi),moduleR*sin(modulePhi),moduleZ);
                 Vector3D nModuleCenter(moduleR*cos(modulePhi),moduleR*sin(modulePhi),-moduleZ);
                 // the rotation matrix of the module
                 Vector3D moduleLocalY(cos(modulePhi),sin(modulePhi),0.);
                 Vector3D pModuleLocalZ(0.,0.,1.); // take different axis to have the same readout direction
                 Vector3D nModuleLocalZ(0.,0.,-1.); // take different axis to have the same readout direction
                 Vector3D nModuleLocalX = moduleLocalY.cross(nModuleLocalZ);
                 Vector3D pModuleLocalX = moduleLocalY.cross(pModuleLocalZ);
                 // local rotation matrices
                 // create the RotationMatrix - negative side
                 RotationMatrix3D nModuleRotation;
                 nModuleRotation.col(0) = nModuleLocalX;
                 nModuleRotation.col(1) = moduleLocalY;
                 nModuleRotation.col(2) = nModuleLocalZ;
                 // create the RotationMatrix - positive side
                 RotationMatrix3D pModuleRotation;
                 pModuleRotation.col(0) = pModuleLocalX;
                 pModuleRotation.col(1) = moduleLocalY;
                 pModuleRotation.col(2) = pModuleLocalZ;
                 // the transforms for the two modules
                 std::shared_ptr<Transform3D> nModuleTransform(new Transform3D(getTransformFromRotTransl(nModuleRotation,nModuleCenter)));
                 std::shared_ptr<Transform3D> pModuleTransform(new Transform3D(getTransformFromRotTransl(pModuleRotation,pModuleCenter)));
                 // create the modules identifier @TODO Idenfier service 
                 Identifier nModuleIdentifier = Identifier(Identifier::value_type(2*imodule));
                 Identifier pModuleIdentifier = Identifier(Identifier::value_type(2*imodule+1));
                 // create the module 
                 GenericDetectorElement* nmodule = new GenericDetectorElement(nModuleIdentifier, nModuleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                 GenericDetectorElement* pmodule = new GenericDetectorElement(pModuleIdentifier, pModuleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                 
                 // the phi cache
                 nneighbourCachePhi.push_back(nmodule); 
                 pneighbourCachePhi.push_back(pmodule);
                 
                 // memory management - we need a detector store to hold them somewhere @TODO add detector store facility
                 m_posnegModule.push_back(nmodule);
                 m_posnegModule.push_back(pmodule);

                 // and the backside one (if configured to do so)
                 if (m_posnegModuleBacksideGap.size()){
                     // ncrease the counter @TODO switch to identifier service
                     nModuleIdentifier = Identifier(Identifier::value_type(++imodule));
                     pModuleIdentifier = Identifier(Identifier::value_type(++imodule));
                     // the new centers
                     nModuleCenter = nModuleCenter + m_posnegModuleBacksideGap[ipnl][ipnR]*nModuleLocalZ;                  
                     pModuleCenter = pModuleCenter + m_posnegModuleBacksideGap[ipnl][ipnR]*pModuleLocalZ;                  
                     // the new transforms
                     nModuleTransform = std::shared_ptr<Transform3D>(new Transform3D(getTransformFromRotTransl(nModuleRotation,nModuleCenter)));
                     pModuleTransform = std::shared_ptr<Transform3D>(new Transform3D(getTransformFromRotTransl(pModuleRotation,pModuleCenter)));
                     // apply the stereo
                     if (m_posnegModuleBacksideStereo.size()){
                         // twist by the stereo angle
                         double stereoBackSide = m_posnegModuleBacksideStereo[ipnl][ipnR];
                         (*nModuleTransform.get()) *= AngleAxis3D(-stereoBackSide, Vector3D::UnitZ());
                         (*pModuleTransform.get()) *= AngleAxis3D(-stereoBackSide, Vector3D::UnitZ());
                         
                     }
                     // everything is set for the next module
                     GenericDetectorElement* bsnmodule = new GenericDetectorElement(nModuleIdentifier, nModuleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                     GenericDetectorElement* bspmodule = new GenericDetectorElement(pModuleIdentifier, pModuleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                     
                     // the backside phi cache
                     nneighbourCachePhiBackside.push_back(bsnmodule); 
                     pneighbourCachePhiBackside.push_back(bspmodule); 
                     
                     // memory management - we need a detector store to hold them somewhere @TODO add detector store facility
                     m_posnegModule.push_back(bsnmodule);
                     m_posnegModule.push_back(bspmodule);
                 }
                 // create the surface 
                 nsVector.push_back(SurfacePosition(&nmodule->surface(), nModuleCenter));
                 psVector.push_back(SurfacePosition(&pmodule->surface(), pModuleCenter));                    
             } 
             // create the phi binned array
             double phiStep = (maxPhi-minPhi)/(ringModulePositionsPhi.size()-1);
             minPhi -= 0.5*phiStep;
             maxPhi += 0.5*phiStep;
             
             // BinUtilities
             BinUtility* nphiBinUtility = new BinUtility(ringModulePositionsPhi.size(),minPhi,maxPhi,closed,binPhi);
             BinUtility* pphiBinUtility = new BinUtility(ringModulePositionsPhi.size(),minPhi,maxPhi,closed,binPhi);
             nRadialSurfaceArrays.push_back( new BinnedArray1D<const Surface* >(nsVector, nphiBinUtility)  );
             pRadialSurfaceArrays.push_back( new BinnedArray1D<const Surface* >(psVector, pphiBinUtility)  );
             
             // push the phi cache into the overall cache
             nneighbourCache.push_back(nneighbourCachePhi); 
             pneighbourCache.push_back(pneighbourCachePhi);
             if (nneighbourCachePhiBackside.size()) nneighbourCacheBackside.push_back(nneighbourCachePhiBackside);
             if (pneighbourCachePhiBackside.size()) pneighbourCacheBackside.push_back(pneighbourCachePhiBackside);
             
         }
                  
         // create the SurfaceArrays
         SurfaceArray* nSurfaceArray = nullptr;
         SurfaceArray* pSurfaceArray = nullptr;
         // we only have a single ring 
         if (nRadialSurfaceArrays.size() == 1 && pRadialSurfaceArrays.size() == 1){
             // just take the one you have
             nSurfaceArray = nRadialSurfaceArrays[0];
             pSurfaceArray = pRadialSurfaceArrays[0];             
         } else {
             // r boundaries ----------------------------------------------
             std::vector< float > rBoundaries = { (float)layerRmin, (float)layerRmax };
             // 
             std::vector< std::pair< SurfaceArray*, Vector3D > > pSurfaceArraysPosition;
             std::vector< std::pair< SurfaceArray*, Vector3D > > nSurfaceArraysPosition;
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
                 pSurfaceArraysPosition.push_back( std::pair< SurfaceArray*, Vector3D>(pRadialSurfaceArrays[irb], Vector3D(0.,0.,innerR)) );
                 nSurfaceArraysPosition.push_back( std::pair< SurfaceArray*, Vector3D>(nRadialSurfaceArrays[irb], Vector3D(0.,0.,innerR)) );
             }
             // and the last one
             pSurfaceArraysPosition.push_back( std::pair< SurfaceArray*, Vector3D>(pRadialSurfaceArrays[layerModuleRadii.size()-1], Vector3D(0.,0.,outerR)) );
             nSurfaceArraysPosition.push_back( std::pair< SurfaceArray*, Vector3D>(nRadialSurfaceArrays[layerModuleRadii.size()-1], Vector3D(0.,0.,outerR)) );
             // sort the rBoundaries for the steering bin 
             std::sort(rBoundaries.begin(), rBoundaries.end());
             // building the 1D 1D 
             BinUtility* nrBinUtility = new BinUtility(rBoundaries, open, binR);
             BinUtility* prBinUtility = new BinUtility(rBoundaries, open, binR);
             // now create the surface arrays
             nSurfaceArray = new BinnedArrayArray< const Surface* >(nSurfaceArraysPosition ,nrBinUtility);
             pSurfaceArray = new BinnedArrayArray< const Surface* >(pSurfaceArraysPosition ,prBinUtility);                         
         }
         
         // create the share disc bounds
         std::shared_ptr<const DiscBounds> dBounds(new RadialBounds(layerRmin-layerEnvelopeR,layerRmax+layerEnvelopeR));
         
         
         MSG_VERBOSE("Filling the neighbour cache for the overlap descriptor.");
         
         // --- interlude --- fill the neighbour cache
         // number of R rings
         size_t rRings = nneighbourCache.size();
         for (size_t iR = 0; iR < rRings; ++iR){
             // loop over the phi values
             size_t phiModules = nneighbourCache[iR].size();
             for (size_t iphi = 0; iphi < phiModules; ++iphi){
                 // create the neighbours
                 std::vector<const DetectorElementBase*> nneighbours;
                 std::vector<const DetectorElementBase*> pneighbours;
                 // closed for phi
                 size_t nphi = (iphi < (phiModules-1)) ? iphi+1 : 0;
                 size_t pphi = (iphi > 0) ? iphi-1 : phiModules-1;
                 // fill neighbours
                 // get the current z with next in phi
                 nneighbours.push_back(nneighbourCache[iR][pphi]);
                 nneighbours.push_back(nneighbourCache[iR][nphi]);
                 pneighbours.push_back(pneighbourCache[iR][pphi]);
                 pneighbours.push_back(pneighbourCache[iR][nphi]);
                 // check if we have backside  modules
                 if (nneighbourCacheBackside.size()){
                     // 
                     nneighbours.push_back(nneighbourCacheBackside[iR][pphi]);
                     nneighbours.push_back(nneighbourCacheBackside[iR][nphi]);
                     pneighbours.push_back(pneighbourCacheBackside[iR][pphi]);
                     pneighbours.push_back(pneighbourCacheBackside[iR][nphi]);                     
                 }
                 
                 // set it for the moment, is missing the ring overlap
                 if (nneighbourCache[iR][iphi]) nneighbourCache[iR][iphi]->registerNeighbours(nneighbours);
                 if (pneighbourCache[iR][iphi]) pneighbourCache[iR][iphi]->registerNeighbours(pneighbours);
                 
                 // if we have the backside module, we should deal with them
                 if (nneighbourCacheBackside.size()){
                     // create the neighbours
                     std::vector<const DetectorElementBase*> bsnneighbours; bsnneighbours.reserve(nneighbours.size());
                     std::vector<const DetectorElementBase*> bspneighbours; bspneighbours.reserve(pneighbours.size());
                     // 
                     for (auto& nnmodule : nneighbours){
                         if (nnmodule != nneighbourCacheBackside[iR][iphi])
                             bsnneighbours.push_back(nnmodule);
                     }
                     for (auto& pnmodule : pneighbours){
                         if (pnmodule != pneighbourCacheBackside[iR][iphi])
                             bspneighbours.push_back(pnmodule);
                     }
                     // now register the final one
                     bsnneighbours.push_back(nneighbourCache[iR][iphi]);
                     bspneighbours.push_back(nneighbourCache[iR][iphi]);
                     nneighbourCacheBackside[iR][iphi]->registerNeighbours(bsnneighbours);
                     pneighbourCacheBackside[iR][iphi]->registerNeighbours(bspneighbours);
                 }
             }
         }
                  
         // layer thickness
         double layerThickness = layerZmax-layerZmin;
         
         // create the layer transforms
         Transform3D* nLayerTransform = new Transform3D(Transform3D::Identity());
         nLayerTransform->translation() = Vector3D(0.,0.,-layerPosZ);
         Transform3D* pLayerTransform = new Transform3D(Transform3D::Identity());
         pLayerTransform->translation() = Vector3D(0.,0.,layerPosZ);

         // create the layers
         LayerPtr nLayer = DiscLayer::create(std::shared_ptr<Transform3D>(nLayerTransform), 
                                             dBounds,
                                             nSurfaceArray,
                                             layerThickness,
                                             nullptr,
                                             nullptr,
                                             active);
                                             
         LayerPtr pLayer = DiscLayer::create(std::shared_ptr<Transform3D>(pLayerTransform), 
                                             dBounds,
                                             pSurfaceArray,
                                             layerThickness,
                                             nullptr,
                                             nullptr,
                                             active);

         // the layer is built le't see if it needs material
         if (m_posnegLayerMaterialProperties.size()){
             // get the material from configuration
             double lMaterialThickness = m_posnegLayerMaterialProperties[ipnl][0];
             double lMaterialX0        = m_posnegLayerMaterialProperties[ipnl][1];
             double lMaterialL0        = m_posnegLayerMaterialProperties[ipnl][2];
             double lMaterialA         = m_posnegLayerMaterialProperties[ipnl][3];
             double lMaterialZ         = m_posnegLayerMaterialProperties[ipnl][4];
             double lMaterialRho       = m_posnegLayerMaterialProperties[ipnl][5];
             MaterialProperties layerMaterialProperties(lMaterialThickness,lMaterialX0,lMaterialL0,lMaterialA,lMaterialZ,lMaterialRho);
             std::shared_ptr<const SurfaceMaterial> layerMaterialPtr(new HomogeneousSurfaceMaterial(layerMaterialProperties));   
             
             // central material
             if (m_posnegLayerMaterialConcentration[ipnl] == 0.){
                 // assign the surface material - the layer surface is the material surface
                 nLayer->surfaceRepresentation().setSurfaceMaterial(layerMaterialPtr);
                 pLayer->surfaceRepresentation().setSurfaceMaterial(layerMaterialPtr);
             } else {
                 // approach surface material
                 // get the approach descriptor - at this stage we know that the approachDescriptor exists
                 auto nApproachSurfaces = nLayer->approachDescriptor()->containedSurfaces();
                 auto pApproachSurfaces = pLayer->approachDescriptor()->containedSurfaces();
                 if (m_posnegLayerMaterialConcentration[ipnl] > 0.){
                     nApproachSurfaces[0]->setSurfaceMaterial(layerMaterialPtr);
                     pApproachSurfaces[1]->setSurfaceMaterial(layerMaterialPtr);
                 } else {
                     nApproachSurfaces[1]->setSurfaceMaterial(layerMaterialPtr);
                     pApproachSurfaces[0]->setSurfaceMaterial(layerMaterialPtr);
                 }
             }
         }                                             
                                             
         // push it into the layer vector
         m_nLayers->push_back(nLayer);
         m_pLayers->push_back(pLayer);
     }
    }
    
    // everything was successful - let's return back
    return StatusCode::SUCCESS;
}

void Acts::GenericLayerBuilder::moduleExtend(const Acts::GenericDetectorElement& detElement,
                                            double thickness, 
                                            double minHalfX, double maxHalfX, double halfY,
                                            double& rmin, double& rmax, 
                                            double& zmin, double& zmax){
    
    // brute force method to find the radial extends for the layer construction
    Vector3D localX = detElement.transform().rotation().col(0);
    Vector3D localY = detElement.transform().rotation().col(1);
    Vector3D localZ = detElement.transform().rotation().col(2);
    // center
    Vector3D center = detElement.transform().translation();    
    // construct the edges 
    Vector3D edge000 = center-0.5*thickness*localZ-minHalfX*localX-halfY*localY;
    Vector3D edge010 = center-0.5*thickness*localZ+minHalfX*localX-halfY*localY;
    Vector3D edge001 = center-0.5*thickness*localZ-maxHalfX*localX+halfY*localY;
    Vector3D edge011 = center-0.5*thickness*localZ+maxHalfX*localX+halfY*localY;
    Vector3D edge100 = center+0.5*thickness*localZ-minHalfX*localX-halfY*localY;
    Vector3D edge110 = center+0.5*thickness*localZ+minHalfX*localX-halfY*localY;
    Vector3D edge101 = center+0.5*thickness*localZ-maxHalfX*localX+halfY*localY;
    Vector3D edge111 = center+0.5*thickness*localZ+maxHalfX*localX+halfY*localY;
    // layer r min / max 
    takeSmallerBigger(rmin,rmax,edge000.perp());
    takeSmallerBigger(rmin,rmax,edge001.perp());
    takeSmallerBigger(rmin,rmax,edge010.perp());
    takeSmallerBigger(rmin,rmax,edge011.perp());
    takeSmallerBigger(rmin,rmax,edge100.perp());
    takeSmallerBigger(rmin,rmax,edge101.perp());
    takeSmallerBigger(rmin,rmax,edge110.perp());
    takeSmallerBigger(rmin,rmax,edge111.perp());
    // layer z min / max
    takeSmallerBigger(zmin,zmax,edge000.z());
    takeSmallerBigger(zmin,zmax,edge001.z());
    takeSmallerBigger(zmin,zmax,edge010.z());
    takeSmallerBigger(zmin,zmax,edge011.z());
    takeSmallerBigger(zmin,zmax,edge100.z());
    takeSmallerBigger(zmin,zmax,edge101.z());
    takeSmallerBigger(zmin,zmax,edge110.z());
    takeSmallerBigger(zmin,zmax,edge111.z());
}


