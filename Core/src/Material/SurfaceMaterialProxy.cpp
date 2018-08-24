// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialProxy.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/SurfaceMaterialProxy.hpp"

Acts::SurfaceMaterialProxy::SurfaceMaterialProxy()
  : SurfaceMaterial(), m_binUtility(nullptr)
{
}

Acts::SurfaceMaterialProxy::SurfaceMaterialProxy(const BinUtility& binUtility)
  : SurfaceMaterial()
  , m_binUtility(std::make_unique<const BinUtility>(binUtility))
{
}

Acts::SurfaceMaterialProxy::SurfaceMaterialProxy(
    const SurfaceMaterialProxy& smproxy)
  : SurfaceMaterial(), m_binUtility(nullptr)
{
  if (smproxy.m_binUtility) {
    m_binUtility = std::make_unique<const BinUtility>(*smproxy.m_binUtility);
  }
}

Acts::SurfaceMaterialProxy*
Acts::SurfaceMaterialProxy::clone() const
{
  return (new SurfaceMaterialProxy(*this));
}

Acts::SurfaceMaterialProxy&
Acts::SurfaceMaterialProxy::operator*=(double)
{
  return (*this);
}

std::ostream&
Acts::SurfaceMaterialProxy::dump(std::ostream& sl) const
{
  sl << "Acts::SurfaceMaterialProxy : " << std::endl;
  if (m_binUtility) {
    sl << "   - Number of Material bins [0,1] : " << m_binUtility->bins(0)
       << " / " << m_binUtility->bins(1) << std::endl;
  } else {
    sl << "   - Homogeneous Material" << std::endl;
  }
  return sl;
}
