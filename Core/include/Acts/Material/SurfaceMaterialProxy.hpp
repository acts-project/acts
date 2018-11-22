// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialProxy.hpp, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Utilities/BinUtility.hpp"

namespace Acts {

/// @class SurfaceMaterialProxy
///
/// @brief proxy to SurfaceMaterial hand over BinUtility
///
/// The SurfaceMaterialProxy class acts as a proxy to the SurfaceMaterial
/// to mark the layers and surfaces on which the material should be mapped on
/// at construction time of the geometry and to hand over the granularity of
/// of the material map with the bin Utility.

class SurfaceMaterialProxy : public SurfaceMaterial
{
public:
  /// Constructor without BinUtility - homogenous material
  SurfaceMaterialProxy() = default;

  /// Constructor with BinUtility - multidimensional material
  ///
  /// @param binUtility a BinUtility determining the granularity
  ///        and binning of the material on the surface/layer
  SurfaceMaterialProxy(const BinUtility& binUtility);

  /// Copy constuctor
  ///
  /// @param smproxy The source proxy
  SurfaceMaterialProxy(const SurfaceMaterialProxy& smproxy) = default;

  /// Copy move constuctor
  ///
  /// @param smproxy The source proxy
  SurfaceMaterialProxy(SurfaceMaterialProxy&& smproxy) = default;

  /// Destructor
  ///
  /// @param smproxy The source proxy
  ~SurfaceMaterialProxy() override = default;

  /// Assignment operator
  ///
  /// @param smproxy The source proxy
  SurfaceMaterialProxy&
  operator=(const SurfaceMaterialProxy& smproxy)
      = default;

  /// Assigment move operator
  ///
  /// @param smproxy The source proxy
  SurfaceMaterialProxy&
  operator=(SurfaceMaterialProxy&& smproxy)
      = default;

  /// Scale operator
  ///
  /// @param
  SurfaceMaterialProxy&
  operator*=(double scale) final;

  /// Return the BinUtility
  const BinUtility&
  binUtility() const;

  /// Return method for full material description of the Surface - from local
  /// coordinates
  ///
  /// @param lp is local positioning vector
  ///
  /// @return will return dummy material
  const MaterialProperties&
  materialProperties(const Vector2D& lp) const final;

  /// Return method for full material description of the Surface - from the
  /// global coordinates
  ///
  /// @param gp is the global positioning vector
  ///
  /// @return will return dummy material
  const MaterialProperties&
  materialProperties(const Vector3D& gp) const final;

  /// Direct access via bins to the MaterialProperties
  ///
  /// @param ib0 indicates the first bin
  /// @param ib1 indicates the seconf bin
  ///
  /// @return will return dummy material
  const MaterialProperties&
  materialProperties(size_t ib0, size_t ib1) const final;

  /// Output Method for std::ostream, to be overloaded by child classes
  std::ostream&
  dump(std::ostream& sl) const final;

private:
  /// two dimensional BinUtility determining
  /// the granularity and binning of the
  /// material on the surface/layer
  BinUtility m_binUtility;

  /// Dummy material properties
  MaterialProperties m_materialProperties;
};
}

inline const Acts::MaterialProperties&
Acts::SurfaceMaterialProxy::materialProperties(const Vector2D& /*lp*/) const
{
  return (m_materialProperties);
}

inline const Acts::MaterialProperties&
Acts::SurfaceMaterialProxy::materialProperties(const Vector3D& /*gp*/) const
{
  return (m_materialProperties);
}

inline const Acts::MaterialProperties&
    Acts::SurfaceMaterialProxy::materialProperties(size_t /*ib0*/,
                                                   size_t /*ib1*/) const
{
  return (m_materialProperties);
}

inline const Acts::BinUtility&
Acts::SurfaceMaterialProxy::binUtility() const
{
  return m_binUtility;
}