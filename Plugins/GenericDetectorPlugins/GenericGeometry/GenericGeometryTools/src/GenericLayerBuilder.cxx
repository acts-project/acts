///////////////////////////////////////////////////////////////////
// GenericLayerBuilder.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core module
#include "Algebra/AlgebraHelper.h"
// Geometry module
#include "GeometryInterfaces/ILayerCreator.h"
#include "GeometryUtils/BinUtility.h"
#include "GeometryUtils/BinnedArray2D.h"
#include "GeometryUtils/BinnedArray1D.h"
#include "GeometryUtils/BinnedArrayArray.h"
#include "GeometryUtils/ApproachDescriptor.h"
#include "DetectorElementBase/DetectorElementBase.h"
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
  m_layerCreator(""),
  m_nLayers(nullptr),
  m_cLayers(nullptr),    
  m_pLayers(nullptr),
  m_approachSurfaceEnvelope(0.5),
  m_centralLayerBinPhimultiplier(1),
  m_centralLayerBinZmultiplier(1),
  m_centralPassiveLayerBuilder(""),
  m_posnegLayerBinRmultiplier(1),  
  m_posnegLayerBinPhimultiplier(1),
  m_posnegPassiveLayerBuilder("")
{
    declareInterface<ILayerBuilder>(this);
   
    //  layer identificaiton
    declareProperty("LayerIdentification",                m_layerIdentification);
    
    declareProperty("LayerCreator",                       m_layerCreator);
   
    // approach surface envelope
    declareProperty("ApproachSurfaceEnvelope",            m_approachSurfaceEnvelope);
   
    // the central layers 
    declareProperty("CentralLayerRadii",                  m_centralLayerRadii);
    declareProperty("CentralLayerEnvelopeR",              m_centralLayerEnvelopeR);
    declareProperty("CentralLayerEnvelopeZ",              m_centralLayerEnvelopeZ);
    declareProperty("CentralLayerModulePhiBinMultiplier", m_centralLayerBinPhimultiplier);
    declareProperty("CentralLayerModuleZbinMultiplier",   m_centralLayerBinZmultiplier);
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
    declareProperty("PosNegLayerModuleRbinMultiplier",    m_posnegLayerBinRmultiplier);
    declareProperty("PosNegLayerModulePhiBinMultiplier",  m_posnegLayerBinPhimultiplier);
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
    
    // retreive the LayerCreator - without that, no layer
    RETRIEVE_FATAL(m_layerCreator);
    
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
           // some screen output
           MSG_VERBOSE("Build layer " << icl << " with target radius = " << layerR);
           // surface vector 
           std::vector<const Surface*> sVector;
           // z/phi values for this layer
           std::vector<double> zValues   = m_centralModulePositionZ[icl];
           std::vector<double> phiValues = m_centralModulePositionPhi[icl];
           std::sort(phiValues.begin(),phiValues.end());
           std::sort(zValues.begin(),zValues.end());
           // envelope stagger & cover
           double layerModuleStaggerZ = m_centralModuleStaggerZ.size() ? m_centralModuleStaggerZ[icl] : 0.;
           double layerEnvelopeCoverZ = m_centralLayerEnvelopeZ.size() ? m_centralLayerEnvelopeZ[icl] : 0.;
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
               // additional z stagger
               double staggerFromZ = stagger ? -0.5*layerModuleStaggerZ : 0.5 * layerModuleStaggerZ; stagger = !stagger;
               // loop of phi values
               for (auto& modulePhi : phiValues){
                   // stagger the modules
                   double moduleR = layerR + staggerFromZ;
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
                   DetectorElementBase* module = new GenericDetectorElement(moduleIdentifier, moduleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                   // register the surface 
                   sVector.push_back(&module->surface());
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
                       DetectorElementBase* bsmodule = new GenericDetectorElement(moduleIdentifier, moduleTransform, moduleBounds, moduleThickness, moduleMaterialPtr);
                       // register the backside as bin member
                       std::vector<const DetectorElementBase*> bsbinmember = {module};
                       std::vector<const DetectorElementBase*> binmember  = {bsmodule};
                       bsmodule->registerBinmembers(bsbinmember);   
                       module->registerBinmembers(binmember);   
                       // memory management - we need a detector store to hold them somewhere @TODO detector store facility
                       m_centralModule.push_back(bsmodule);
                   }
               }
           }
           
           size_t phiBins = m_centralLayerBinPhimultiplier*phiValues.size();
           size_t zBins   = m_centralLayerBinZmultiplier*zValues.size();
           
           // create the surface array - it will also fill the accesible binmember chache if avalable
           LayerPtr cLayer = m_layerCreator->cylinderLayer(sVector, m_approachSurfaceEnvelope, layerEnvelopeCoverZ, phiBins, zBins);               
                      
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
         size_t layerBinsR                         = 0;
         size_t layerBinsPhi                       = 0;
             
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
         std::vector< const Surface* > nsVector;
         std::vector< const Surface* > psVector;
         
         // loop over bins in R
         size_t imodule = 0;

         // staggering sterring
         bool rstagger = true;
         // screen output
         MSG_VERBOSE("This pair of discs has " << layerModuleRadii.size() << " rings.");
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
             if (layerModuleMinHalfX.size() && moduleMinHalfX != moduleMaxHalfX)
                   pBounds = new TrapezoidBounds(moduleMinHalfX, moduleMaxHalfX, moduleHalfY);
             else 
                   pBounds = new RectangleBounds(moduleMaxHalfX, moduleHalfY); 
             // now create the shared bounds from it
             std::shared_ptr<const PlanarBounds> moduleBounds(pBounds);
             // stagger in phi
             bool phistagger = true;
             double modulePhiStagger = layerModulePhiStagger[ipnR];

             // the phi module of this ring
             auto ringModulePositionsPhi = layerModulePositionsPhi[ipnR];
             
             MSG_VERBOSE("Ring - " << ipnR << " - has " << ringModulePositionsPhi.size() << " phi modules.");
             
             // now loop over phi
             for (auto& modulePhi : ringModulePositionsPhi){
                 // update the module z position
                 moduleZ += phistagger ? 0.5*modulePhiStagger : -0.5*modulePhiStagger; phistagger = !phistagger;
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
                     
                     // register the backside of the binmembers
                     std::vector<const DetectorElementBase*> bspbinmember = {pmodule};
                     std::vector<const DetectorElementBase*> pbinmember = {bspmodule};
                     std::vector<const DetectorElementBase*> bsnbinmember = {nmodule};
                     std::vector<const DetectorElementBase*> nbinmember = {bsnmodule};
                     bsnmodule->registerBinmembers(bsnbinmember);   
                     nmodule->registerBinmembers(nbinmember);   
                     bspmodule->registerBinmembers(bspbinmember);   
                     pmodule->registerBinmembers(pbinmember);   
                     
                     // memory management - we need a detector store to hold them somewhere @TODO add detector store facility
                     m_posnegModule.push_back(bsnmodule);
                     m_posnegModule.push_back(bspmodule);
                 }
                 // create the surface 
                 nsVector.push_back(&nmodule->surface());
                 psVector.push_back(&pmodule->surface());                    
             } 
             // take the values of the maximum phi bins
             if (ringModulePositionsPhi.size() > layerBinsPhi){
                 layerBinsPhi = ringModulePositionsPhi.size();
             }
                          
         }
         
         // estimate teh layerBinsR 
         layerBinsR = layerModuleRadii.size() > 1 ? layerModuleRadii.size()*m_posnegLayerBinRmultiplier : 1;
         layerBinsPhi *= m_posnegLayerBinPhimultiplier;
         
         // create teh surface arrays 
         LayerPtr nLayer = m_layerCreator->discLayer(nsVector, layerEnvelopeR, layerEnvelopeR, m_approachSurfaceEnvelope, layerBinsR, layerBinsPhi);
         LayerPtr pLayer = m_layerCreator->discLayer(psVector, layerEnvelopeR, layerEnvelopeR, m_approachSurfaceEnvelope, layerBinsR, layerBinsPhi);

         // create the layer transforms
         Transform3D* nLayerTransform = new Transform3D(Transform3D::Identity());
         nLayerTransform->translation() = Vector3D(0.,0.,-layerPosZ);
         Transform3D* pLayerTransform = new Transform3D(Transform3D::Identity());
         pLayerTransform->translation() = Vector3D(0.,0.,layerPosZ);

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
