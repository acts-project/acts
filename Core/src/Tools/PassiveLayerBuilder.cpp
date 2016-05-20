///////////////////////////////////////////////////////////////////
// PassiveLayerBuilder.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/PassiveLayerBuilder.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Material/MaterialProperties.hpp"

// constructor
Acts::PassiveLayerBuilder::PassiveLayerBuilder(const PassiveLayerBuilder::Config& plConfig) :
  m_config(),
  m_nLayers(),
  m_cLayers(),
  m_pLayers(),
  m_constructionFlag(false)
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
        ACTS_DEBUG("Configured to build " << numcLayers << " passive central layers.");
        m_cLayers.reserve(numcLayers);
        // loop through
        for (size_t icl = 0; icl < numcLayers; ++icl){
            // some screen output
            ACTS_VERBOSE("- build layer " << icl << " with radius = " << m_config.centralLayerRadii.at(icl) << " and halfZ = " << m_config.centralLayerHalflengthZ.at(icl));
            // create the layer and push it back
            auto cBounds = std::make_shared<CylinderBounds>(m_config.centralLayerRadii[icl],m_config.centralLayerHalflengthZ.at(icl));
            // create the layer
            LayerPtr cLayer = CylinderLayer::create(nullptr, cBounds, nullptr, m_config.centralLayerThickness.at(icl));
            // assign the material to the layer surface
            std::shared_ptr<const SurfaceMaterial> material = nullptr;
            // create the material from jobOptions
            if (m_config.centralLayerMaterialX0.size()){
                // create homogeneous material
                material = std::make_shared<const HomogeneousSurfaceMaterial>(MaterialProperties(m_config.centralLayerThickness.at(icl),
                                                                                                 m_config.centralLayerMaterialX0.at(icl),
                                                                                                 m_config.centralLayerMaterialL0.at(icl),
                                                                                                 m_config.centralLayerMaterialA.at(icl),
                                                                                                 m_config.centralLayerMaterialZ.at(icl),
                                                                                                 m_config.centralLayerMaterialRho.at(icl)), 1.);

                // sign it to the surface
                cLayer->surfaceRepresentation().setSurfaceMaterial(material);
            }
            // push it into the layer vector
            m_cLayers.push_back(cLayer);
        }
    }

    // pos/neg layers
    size_t numpnLayers = m_config.posnegLayerPositionZ.size();
    if (numpnLayers){
        ACTS_DEBUG("Configured to build 2 * " << numpnLayers << " passive positive/negative side layers.");
        m_pLayers.reserve(numpnLayers);
        m_nLayers.reserve(numpnLayers);
        // loop through
        for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl){
            // some screen output
            ACTS_VERBOSE("- build layers " << (2*ipnl) << " and "<<  (2*ipnl)+1 << " at +/- z = " << m_config.posnegLayerPositionZ.at(ipnl));
            //                               << " and rMin/rMax = " << m_config.posnegLayerRmin.at(ipnl) << " / " << m_config.posnegLayerRmax.at(ipnl));
            // create the share disc bounds
            std::shared_ptr<const DiscBounds> dBounds = std::make_shared<RadialBounds>(m_config.posnegLayerRmin.at(ipnl), m_config.posnegLayerRmax.at(ipnl));
            // create the layer transforms
            Transform3D* nTransform = new Transform3D(Transform3D::Identity());
            nTransform->translation() = Vector3D(0.,0.,-m_config.posnegLayerPositionZ.at(ipnl));
            Transform3D* pTransform = new Transform3D(Transform3D::Identity());
            pTransform->translation() = Vector3D(0.,0.,m_config.posnegLayerPositionZ.at(ipnl));
            // create the layers
            LayerPtr nLayer = DiscLayer::create(std::shared_ptr<Transform3D>(nTransform), dBounds, nullptr, m_config.posnegLayerThickness.at(ipnl));
            LayerPtr pLayer = DiscLayer::create(std::shared_ptr<Transform3D>(pTransform), dBounds, nullptr, m_config.posnegLayerThickness.at(ipnl));
            // assign the material to the layer surface
            std::shared_ptr<const SurfaceMaterial> material = nullptr;
            // create the material from jobOptions
            if (m_config.posnegLayerMaterialX0.size()){
                // create homogeneous material
                material = std::make_shared<const HomogeneousSurfaceMaterial>(MaterialProperties(m_config.posnegLayerThickness.at(ipnl),
                                                                                                                    m_config.posnegLayerMaterialX0.at(ipnl),
                                                                                                                    m_config.posnegLayerMaterialL0.at(ipnl),
                                                                                                                    m_config.posnegLayerMaterialA.at(ipnl),
                                                                                                                    m_config.posnegLayerMaterialZ.at(ipnl),
                                                                                                 m_config.posnegLayerMaterialRho.at(ipnl)),1.);
                // sign it to the surface
                nLayer->surfaceRepresentation().setSurfaceMaterial(material);
                pLayer->surfaceRepresentation().setSurfaceMaterial(material);
            }
            // push it into the layer vector
            m_nLayers.push_back(nLayer);
            m_pLayers.push_back(pLayer);
        }
    }

    m_constructionFlag = true;
    return true;
}
