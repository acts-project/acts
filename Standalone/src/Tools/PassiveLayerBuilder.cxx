///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core module
#include "Algebra/AlgebraDefinitions.h"
// Geometry module
#include "GeometryTools/PassiveLayerBuilder.h"
#include "Detector/CylinderLayer.h"
#include "Detector/DiscLayer.h"
#include "Surfaces/CylinderBounds.h"
#include "Surfaces/RadialBounds.h"
#include "Material/HomogeneousSurfaceMaterial.h"
#include "Material/MaterialProperties.h"


DECLARE_TOOL_FACTORY(Acts::PassiveLayerBuilder)

// constructor
Acts::PassiveLayerBuilder::PassiveLayerBuilder(const std::string& t, const std::string& n, const IInterface* p) :
  Acts::AlgToolBase(t,n,p),
  m_layerIdentification("PassiveLayerBuilder"),   
  m_nLayers(nullptr),
  m_cLayers(nullptr),    
  m_pLayers(nullptr)
{
    declareInterface<ILayerBuilder>(this);   
    // read in the layer specs from jobOptions
    // the central layers 
    declareProperty("CentralLayerRadii",          m_centralLayerRadii);
    declareProperty("CentralLayerHalflengthZ",    m_centralLayerHalflengthZ);
    declareProperty("CentralLayerThickness",      m_centralLayerThickness);
    declareProperty("CentralLayerMaterialX0",     m_centralLayerMaterialX0);
    declareProperty("CentralLayerMaterialL0",     m_centralLayerMaterialL0);
    declareProperty("CentralLayerMaterialA",      m_centralLayerMaterialA);
    declareProperty("CentralLayerMaterialZ",      m_centralLayerMaterialZ);
    declareProperty("CentralLayerMaterialRho",    m_centralLayerMaterialRho);
    // the layers at p/n side 
    declareProperty("PosnegLayerPositionZ",       m_posnegLayerPositionZ);
    declareProperty("PosnegLayerRmin",            m_posnegLayerRmin);
    declareProperty("PosnegLayerRmax",            m_posnegLayerRmax);
    declareProperty("PosnegLayerThickness",       m_posnegLayerThickness);
    declareProperty("PosnegLayerMaterialX0",      m_posnegLayerMaterialX0);
    declareProperty("PosnegLayerMaterialL0",      m_posnegLayerMaterialL0);
    declareProperty("PosnegLayerMaterialA",       m_posnegLayerMaterialA);
    declareProperty("PosnegLayerMaterialZ",       m_posnegLayerMaterialZ);
    declareProperty("PosnegLayerMaterialRho",     m_posnegLayerMaterialRho);
    // layer identification given by job options
    declareProperty("LayerIdentification",        m_layerIdentification);
}

// destructor
Acts::PassiveLayerBuilder::~PassiveLayerBuilder()
{}

// initialize
StatusCode Acts::PassiveLayerBuilder::initialize()
{
    MSG_DEBUG( "initialize()" );
    //Tool needs to be initialized
    if (!AlgToolBase::initialize()) return StatusCode::FAILURE;
    return constructLayers();
}

//finalize
StatusCode Acts::PassiveLayerBuilder::finalize()
{
    MSG_DEBUG( "finalize()" );
    return StatusCode::SUCCESS;
}


StatusCode Acts::PassiveLayerBuilder::constructLayers() 
{
    
    // the central layers
    size_t numcLayers = m_centralLayerRadii.size();
    if (numcLayers){
        MSG_DEBUG("Configured to build " << numcLayers << " passive central layers.");
        m_cLayers = new Acts::LayerVector;
        m_cLayers->reserve(numcLayers);
        // loop through
        for (size_t icl = 0; icl < numcLayers; ++icl){
            // some screen output
            MSG_VERBOSE("- build layer " << icl << " with radius = " << m_centralLayerRadii[icl] << " and halfZ = " << m_centralLayerHalflengthZ[icl]);
            // create the layer and push it back
            std::shared_ptr<const CylinderBounds> cBounds(new CylinderBounds(m_centralLayerRadii[icl],m_centralLayerHalflengthZ[icl]));
            // create the layer
            LayerPtr cLayer = CylinderLayer::create(nullptr, cBounds, nullptr, m_centralLayerThickness[icl]);
            // assign the material to the layer surface
            std::shared_ptr<const SurfaceMaterial> material = nullptr;
            // create the material from jobOptions
            if (m_centralLayerMaterialX0.size()){
                // create homogeneous material
                material = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(MaterialProperties(m_centralLayerThickness[icl],
                                                                                                                    m_centralLayerMaterialX0[icl],
                                                                                                                    m_centralLayerMaterialL0[icl],
                                                                                                                    m_centralLayerMaterialA[icl],
                                                                                                                    m_centralLayerMaterialZ[icl],
                                                                                                                    m_centralLayerMaterialRho[icl]), 1.));
                // sign it to the surface
                cLayer->surfaceRepresentation().setSurfaceMaterial(material);
            } 
            // push it into the layer vector
            m_cLayers->push_back(cLayer);
        }
    }
    
    // pos/neg layers
    size_t numpnLayers = m_posnegLayerPositionZ.size();
    if (numpnLayers){
        MSG_DEBUG("Configured to build 2 * " << numpnLayers << " passive positive/negative side layers.");
        m_pLayers = new Acts::LayerVector;
        m_pLayers->reserve(numpnLayers);
        m_nLayers = new Acts::LayerVector;
        m_nLayers->reserve(numpnLayers);
        // loop through
        for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl){
            // some screen output
            MSG_VERBOSE("- build layers " << (2*ipnl) << " and "<<  (2*ipnl)+1 << " at +/- z = " << m_posnegLayerPositionZ[ipnl] 
                                          << " and rMin/rMax = " << m_posnegLayerRmin[ipnl] << " / " << m_posnegLayerRmax[ipnl]);
            // create the share disc bounds
            std::shared_ptr<const DiscBounds> dBounds(new RadialBounds(m_posnegLayerRmin[ipnl], m_posnegLayerRmax[ipnl]));
            // create the layer transforms
            Transform3D* nTransform = new Transform3D(Transform3D::Identity());
            nTransform->translation() = Vector3D(0.,0.,-m_posnegLayerPositionZ[ipnl]);
            Transform3D* pTransform = new Transform3D(Transform3D::Identity());
            pTransform->translation() = Vector3D(0.,0.,m_posnegLayerPositionZ[ipnl]);
            // create the layers
            LayerPtr nLayer = DiscLayer::create(std::shared_ptr<Transform3D>(nTransform), dBounds, nullptr, m_posnegLayerThickness[ipnl]);
            LayerPtr pLayer = DiscLayer::create(std::shared_ptr<Transform3D>(pTransform), dBounds, nullptr, m_posnegLayerThickness[ipnl]);
            // assign the material to the layer surface
            std::shared_ptr<const SurfaceMaterial> material = nullptr;
            // create the material from jobOptions
            if (m_posnegLayerMaterialX0.size()){
                // create homogeneous material
                material = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(MaterialProperties(m_posnegLayerThickness[ipnl],
                                                                                                                    m_posnegLayerMaterialX0[ipnl],
                                                                                                                    m_posnegLayerMaterialL0[ipnl],
                                                                                                                    m_posnegLayerMaterialA[ipnl],
                                                                                                                    m_posnegLayerMaterialZ[ipnl],
                                                                                                                    m_posnegLayerMaterialRho[ipnl]), 1.));
                // sign it to the surface
                nLayer->surfaceRepresentation().setSurfaceMaterial(material);
                pLayer->surfaceRepresentation().setSurfaceMaterial(material);
            } 
            // push it into the layer vector
            m_nLayers->push_back(nLayer);
            m_pLayers->push_back(pLayer);
        }
    }
    return StatusCode::SUCCESS;
}
