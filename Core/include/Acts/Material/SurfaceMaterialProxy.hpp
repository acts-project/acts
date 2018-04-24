// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialProxy.hpp, ACTS project MaterialPlugins
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
  SurfaceMaterialProxy();

  /// Constructor with BinUtility - multidimensional material
  ///
  /// @param binUtility a BinUtility determining the granularity
  ///        and binning of the material on the surface/layer
  SurfaceMaterialProxy(const BinUtility& binUtility);

  /// Copy constuctor
  SurfaceMaterialProxy(const SurfaceMaterialProxy& smproxy);

  /// Destructor
  virtual ~SurfaceMaterialProxy() = default;

  /// Implicit constructor
  /// - uses the copy constructor
  SurfaceMaterialProxy*
  clone() const final override;

  /// Scale operator
  virtual SurfaceMaterialProxy&
  operator*=(double scale) final override;

  /// Return the BinUtility - can be nullptr
  const BinUtility*
  binUtility() const;

  /// Return method for full material description of the Surface - from local
  /// coordinates
  ///
  /// @param lp is local positioning vector
  ///
  /// @return will return dummy material
  virtual const MaterialProperties*
  material(const Vector2D& lp) const final override;

  /// Return method for full material description of the Surface - from the
  /// global coordinates
  ///
  /// @param gp is the global positioning vector
  ///
  /// @return will return dummy material
  virtual const MaterialProperties*
  material(const Vector3D& gp) const final override;

  /// Direct access via bins to the MaterialProperties
  ///
  /// @param ib0 indicates the first bin
  /// @param ib1 indicates the seconf bin
  ///
  /// @return will return dummy material
  virtual const MaterialProperties*
  material(size_t ib0, size_t ib1) const final override;

  /// Output Method for std::ostream, to be overloaded by child classes
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

private:
  /// two dimensional BinUtility determining
  /// the granularity and binning of the
  /// material on the surface/layer
  std::unique_ptr<const BinUtility> m_binUtility;

  /// Dummy material properties
  MaterialProperties m_materialProperties;
};
}

inline const Acts::MaterialProperties*
Acts::SurfaceMaterialProxy::material(const Vector2D&) const
{
  return (&m_materialProperties);
}

inline const Acts::MaterialProperties*
Acts::SurfaceMaterialProxy::material(const Vector3D&) const
{
  return (&m_materialProperties);
}

inline const Acts::MaterialProperties*
    Acts::SurfaceMaterialProxy::material(size_t, size_t) const
{
  return (&m_materialProperties);
}

inline const Acts::BinUtility*
Acts::SurfaceMaterialProxy::binUtility() const
{
  return m_binUtility.get();
}