///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/PassiveLayerBuilder.h"
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Layers/CylinderLayer.h"
#include "ACTS/Layers/DiscLayer.h"
#include "ACTS/Surfaces/CylinderBounds.h"
#include "ACTS/Surfaces/RadialBounds.h"
#include "ACTS/Material/HomogeneousSurfaceMaterial.h"
#include "ACTS/Material/MaterialProperties.h"

// constructor
Acts::PassiveLayerBuilder::PassiveLayerBuilder(const PassiveLayerBuilder::Config& plConfig) :
  m_config(),
  m_constructionFlag(false),
  m_nLayers(nullptr),
  m_cLayers(nullptr),    
  m_pLayers(nullptr)
{
    setConfiguration(plConfig);
}

// configuration method
void Acts::PassiveLayerBuilder::setConfiguration(const PassiveLayerBuilder::Config& plConfig)
{
    //!< @TODO add configuration check
    m_config = plConfig;
}

bool Acts::PassiveLayerBuilder::constructLayers() const
{
    
    // the central layers
    size_t numcLayers = m_config.centralLayerRadii.size();
    if (numcLayers){
        // MSG_DEBUG("Configured to build " << numcLayers << " passive central layers.");
        m_cLayers = new Acts::LayerVector;
        m_cLayers->reserve(numcLayers);
        // loop through
        for (size_t icl = 0; icl < numcLayers; ++icl){
            // some screen output
            // MSG_VERBOSE("- build layer " << icl << " with radius = " << m_config.centralLayerRadii[icl] << " and halfZ = " << m_config.centralLayerHalflengthZ[icl]);
            // create the layer and push it back
            std::shared_ptr<const CylinderBounds> cBounds(new CylinderBounds(m_config.centralLayerRadii[icl],m_config.centralLayerHalflengthZ[icl]));
            // create the layer
            LayerPtr cLayer = CylinderLayer::create(nullptr, cBounds, nullptr, m_config.centralLayerThickness[icl]);
            // assign the material to the layer surface
            std::shared_ptr<const SurfaceMaterial> material = nullptr;
            // create the material from jobOptions
            if (m_config.centralLayerMaterialX0.size()){
                // create homogeneous material
                material = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(MaterialProperties(m_config.centralLayerThickness[icl],
                                                                                                                    m_config.centralLayerMaterialX0[icl],
                                                                                                                    m_config.centralLayerMaterialL0[icl],
                                                                                                                    m_config.centralLayerMaterialA[icl],
                                                                                                                    m_config.centralLayerMaterialZ[icl],
                                                                                                                    m_config.centralLayerMaterialRho[icl]), 1.));
                // sign it to the surface
                cLayer->surfaceRepresentation().setSurfaceMaterial(material);
            } 
            // push it into the layer vector
            m_cLayers->push_back(cLayer);
        }
    }
    
    // pos/neg layers
    size_t numpnLayers = m_config.posnegLayerPositionZ.size();
    if (numpnLayers){
        // MSG_DEBUG("Configured to build 2 * " << numpnLayers << " passive positive/negative side layers.");
        m_pLayers = new Acts::LayerVector;
        m_pLayers->reserve(numpnLayers);
        m_nLayers = new Acts::LayerVector;
        m_nLayers->reserve(numpnLayers);
        // loop through
        for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl){
            // some screen output
            // MSG_VERBOSE("- build layers " << (2*ipnl) << " and "<<  (2*ipnl)+1 << " at +/- z = " << m_config.posnegLayerPositionZ[ipnl] 
            //                               << " and rMin/rMax = " << m_config.posnegLayerRmin[ipnl] << " / " << m_config.posnegLayerRmax[ipnl]);
            // create the share disc bounds
            std::shared_ptr<const DiscBounds> dBounds(new RadialBounds(m_config.posnegLayerRmin[ipnl], m_config.posnegLayerRmax[ipnl]));
            // create the layer transforms
            Transform3D* nTransform = new Transform3D(Transform3D::Identity());
            nTransform->translation() = Vector3D(0.,0.,-m_config.posnegLayerPositionZ[ipnl]);
            Transform3D* pTransform = new Transform3D(Transform3D::Identity());
            pTransform->translation() = Vector3D(0.,0.,m_config.posnegLayerPositionZ[ipnl]);
            // create the layers
            LayerPtr nLayer = DiscLayer::create(std::shared_ptr<Transform3D>(nTransform), dBounds, nullptr, m_config.posnegLayerThickness[ipnl]);
            LayerPtr pLayer = DiscLayer::create(std::shared_ptr<Transform3D>(pTransform), dBounds, nullptr, m_config.posnegLayerThickness[ipnl]);
            // assign the material to the layer surface
            std::shared_ptr<const SurfaceMaterial> material = nullptr;
            // create the material from jobOptions
            if (m_config.posnegLayerMaterialX0.size()){
                // create homogeneous material
                material = std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(MaterialProperties(m_config.posnegLayerThickness[ipnl],
                                                                                                                    m_config.posnegLayerMaterialX0[ipnl],
                                                                                                                    m_config.posnegLayerMaterialL0[ipnl],
                                                                                                                    m_config.posnegLayerMaterialA[ipnl],
                                                                                                                    m_config.posnegLayerMaterialZ[ipnl],
                                                                                                                    m_config.posnegLayerMaterialRho[ipnl]), 1.));
                // sign it to the surface
                nLayer->surfaceRepresentation().setSurfaceMaterial(material);
                pLayer->surfaceRepresentation().setSurfaceMaterial(material);
            } 
            // push it into the layer vector
            m_nLayers->push_back(nLayer);
            m_pLayers->push_back(pLayer);
        }
    }

    m_constructionFlag = true;
    return true;
}
