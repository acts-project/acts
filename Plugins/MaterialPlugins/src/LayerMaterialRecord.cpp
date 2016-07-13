#include "ACTS/Plugins/MaterialPlugins/LayerMaterialRecord.hpp"

Acts::LayerMaterialRecord::LayerMaterialRecord() :
m_binUtility(nullptr),
m_materialMatrix(),
m_entries(0),
m_sumThickness(0.),
m_sumRho(0.)
{}

Acts::LayerMaterialRecord::LayerMaterialRecord(const BinUtility* binutility) :
m_binUtility(binutility),
m_materialMatrix(),
m_entries(0),
m_sumThickness(0.),
m_sumRho(0.)
{}

Acts::LayerMaterialRecord::LayerMaterialRecord(const LayerMaterialRecord& lmrecord) :
m_binUtility(lmrecord.m_binUtility),
m_materialMatrix(lmrecord.m_materialMatrix),
m_entries(lmrecord.m_entries),
m_sumThickness(lmrecord.m_sumThickness),
m_sumRho(lmrecord.m_sumRho)
{}

Acts::LayerMaterialRecord* Acts::LayerMaterialRecord::clone() const
{
    return (new LayerMaterialRecord(*this));
}

Acts::LayerMaterialRecord& Acts::LayerMaterialRecord::operator=(const LayerMaterialRecord& lmrecord)
{
    if (this != &lmrecord) {
        m_binUtility        = lmrecord.m_binUtility;
        m_materialMatrix    = lmrecord.m_materialMatrix;
        m_entries           = lmrecord.m_entries;
        m_sumThickness      = lmrecord.m_sumThickness;
        m_sumRho            = lmrecord.m_sumRho;
    }
    return (*this);
}

void Acts::LayerMaterialRecord::addLayerMaterialProperties(const Acts::Vector3D& pos, const Acts::MaterialProperties* newMaterial)
{
    ++m_entries;
    // get the bins corresponding to the position
    size_t bin1 = m_binUtility->bin(pos,0);
    size_t bin2 = m_binUtility->bin(pos,1);
    
    // get the material which might be there already and add new material and weigh it
    const Acts::MaterialProperties* material = m_materialMatrix.at(bin1).at(bin2);
    float newThickness     = newMaterial->thickness();
    float newRho           = newMaterial->averageRho();
    
    float thickness        = material->thickness() + newThickness;
    m_sumThickness          += newThickness;
    
    float rho              = material->averageRho() + newRho * newThickness;
    m_sumRho                += newRho;
    
    float x0               = material->x0() + newMaterial->x0() * newThickness;
    float l0               = material->l0() + newMaterial->l0() * newThickness;
    float A                = material->averageA() + newMaterial->averageA() * newRho;
    float Z                = material->averageZ() + newMaterial->averageZ() * newRho;
    // set the new current material
    const Acts::Material updatedMaterial(x0,l0,A,Z,rho);
    m_materialMatrix.at(bin1).at(bin2)->setMaterial(updatedMaterial,thickness);
}

void Acts::LayerMaterialRecord::averageMaterial()
{
    if (m_entries > 0) {
        // first finalize the material by averaging over it
        size_t bins1 = m_binUtility->bins(0);
        size_t bins2 = m_binUtility->bins(1);
    
        for (size_t bin1 = 0; bin1 < bins1; bin1++) {
            for (size_t bin2 = 0; bin2 < bins2; bin2++) {
                const Acts::MaterialProperties* material = m_materialMatrix.at(bin1).at(bin2);
                float thickness    = material->thickness()/m_entries;
                float rho          = material->averageRho()/m_sumThickness;
                float x0           = material->x0()/m_sumThickness;
                float l0           = material->l0()/m_sumThickness;
                float A            = material->averageA()/m_sumRho;
                float Z            = material->averageZ()/m_sumRho;
                // set the new current material
                const Acts::Material updatedMaterial(x0,l0,A,Z,rho);
                m_materialMatrix.at(bin1).at(bin2)->setMaterial(updatedMaterial,thickness);
            } // b2
        } // b1
    }
    // reset counter and sums
    m_entries = 0;
    m_sumThickness = 0.;
    m_sumRho = 0.;
}

std::shared_ptr<const Acts::BinnedSurfaceMaterial> Acts::LayerMaterialRecord::layerMaterial()
{
    averageMaterial();
    // return the binned surface material
    return (std::make_shared<const Acts::BinnedSurfaceMaterial>(*m_binUtility,m_materialMatrix));
}
