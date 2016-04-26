///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core module
#include "ACTS/Tools/PassiveLayerBuilder.h"

#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Layers/CylinderLayer.h"
#include "ACTS/Layers/DiscLayer.h"
#include "ACTS/Surfaces/CylinderBounds.h"
#include "ACTS/Surfaces/RadialBounds.h"
#include "ACTS/Material/HomogeneousSurfaceMaterial.h"
#include "ACTS/Material/MaterialProperties.h"


// constructor
Acts::PassiveLayerBuilder::PassiveLayerBuilder() :
  m_constructionFlag(false),
  m_layerIdentification("PassiveLayerBuilder"),   
  m_nLayers(nullptr),
  m_cLayers(nullptr),    
  m_pLayers(nullptr)
{
}

bool Acts::PassiveLayerBuilder::constructLayers() const
{
    
    // the central layers
    size_t numcLayers = m_centralLayerRadii.size();
    if (numcLayers){
        // MSG_DEBUG("Configured to build " << numcLayers << " passive central layers.");
        m_cLayers = new Acts::LayerVector;
        m_cLayers->reserve(numcLayers);
        // loop through
        for (size_t icl = 0; icl < numcLayers; ++icl){
            // some screen output
            // MSG_VERBOSE("- build layer " << icl << " with radius = " << m_centralLayerRadii[icl] << " and halfZ = " << m_centralLayerHalflengthZ[icl]);
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
        // MSG_DEBUG("Configured to build 2 * " << numpnLayers << " passive positive/negative side layers.");
        m_pLayers = new Acts::LayerVector;
        m_pLayers->reserve(numpnLayers);
        m_nLayers = new Acts::LayerVector;
        m_nLayers->reserve(numpnLayers);
        // loop through
        for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl){
            // some screen output
            // MSG_VERBOSE("- build layers " << (2*ipnl) << " and "<<  (2*ipnl)+1 << " at +/- z = " << m_posnegLayerPositionZ[ipnl] 
            //                               << " and rMin/rMax = " << m_posnegLayerRmin[ipnl] << " / " << m_posnegLayerRmax[ipnl]);
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

    m_constructionFlag = true;
    return true;
}
