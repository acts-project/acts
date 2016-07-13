#include "ACTS/Plugins/MaterialPlugins/SurfaceMaterialProxy.hpp"

Acts::SurfaceMaterialProxy::SurfaceMaterialProxy(BinUtility& binutility) :
Acts::SurfaceMaterial(),
m_binUtility(binutility.clone())
{}

Acts::SurfaceMaterialProxy::SurfaceMaterialProxy(const SurfaceMaterialProxy& smproxy) :
Acts::SurfaceMaterial(),
m_binUtility(smproxy.m_binUtility->clone())
{}

Acts::SurfaceMaterialProxy* Acts::SurfaceMaterialProxy::clone() const
{
    return (new SurfaceMaterialProxy(*this));
}

Acts::SurfaceMaterial& Acts::SurfaceMaterialProxy::operator*=(double)
{
    return (*this);
}

std::ostream& Acts::SurfaceMaterialProxy::dump( std::ostream& sl) const
{
    sl << "Acts::SurfaceMaterialProxy : " << std::endl;
    sl << "   - Number of Material bins (1/2) : " << m_binUtility->max(0)+1 << " / " << m_binUtility->max(1)+1  << std::endl;
    return sl;
}