// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialProxy.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Material/SurfaceMaterialProxy.hpp"

Acts::SurfaceMaterialProxy::SurfaceMaterialProxy(BinUtility& binutility)
  : Acts::SurfaceMaterial()
  , m_binUtility(binutility)
{
}

Acts::SurfaceMaterialProxy::SurfaceMaterialProxy(
    const SurfaceMaterialProxy& smproxy)
  : Acts::SurfaceMaterial()
  , m_binUtility(smproxy.m_binUtility)
{
}

Acts::SurfaceMaterialProxy*
Acts::SurfaceMaterialProxy::clone() const
{
  return (new SurfaceMaterialProxy(*this));
}

Acts::SurfaceMaterial&
Acts::SurfaceMaterialProxy::operator*=(double)
{
  return (*this);
}

std::ostream&
Acts::SurfaceMaterialProxy::dump(std::ostream& sl) const
{
  sl << "Acts::SurfaceMaterialProxy : " << std::endl;
  sl << "   - Number of Material bins [0,1] : " << m_binUtility.max(0) + 1
     << " / " << m_binUtility.max(1) + 1 << std::endl;
  return sl;
}