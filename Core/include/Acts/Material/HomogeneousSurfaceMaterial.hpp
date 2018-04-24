// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// HomogeneousSurfaceMaterial.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIAL_HOMOGENOUSLAYERMATERIAL_H
#define ACTS_MATERIAL_HOMOGENOUSLAYERMATERIAL_H

#include <vector>
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class HomogeneousSurfaceMaterial
///
/// It extends the SurfaceMaterial base class and describes a simple homogeneous
/// material descriptions

class HomogeneousSurfaceMaterial : public SurfaceMaterial
{
public:
  /// Default Constructor - deleted
  HomogeneousSurfaceMaterial() = delete;

  /// Explizit constructor
  ///
  /// @param fullmat are the full material properties
  /// @param splitFactor is the split for pre/post update
  HomogeneousSurfaceMaterial(const MaterialProperties& fullmat,
                             double                    splitFactor = 1.);

  /// Copy Constructor
  ///
  /// @param hsm is the source material
  HomogeneousSurfaceMaterial(const HomogeneousSurfaceMaterial& hsm);

  /// Destructor
  virtual ~HomogeneousSurfaceMaterial();

  /// Pseudo-Constructor clone(
  HomogeneousSurfaceMaterial*
  clone() const final override;

  /// Assignment operator
  HomogeneousSurfaceMaterial&
  operator=(const HomogeneousSurfaceMaterial& lmp);

  /// Scale operator
  /// - it is effectively a thickness scaling
  ///
  /// @param scale is the scale factor
  HomogeneousSurfaceMaterial&
  operator*=(double scale) final override;

  /// @copydoc SurfaceMaterial::material(const Vector2D&)
  ///
  /// @note the input parameter is ignored
  virtual const MaterialProperties*
  material(const Vector2D& lp) const final override;

  /// @copydoc SurfaceMaterial::material(const Vector3D&)
  ///
  /// @note the input parameter is ignored
  virtual const MaterialProperties*
  material(const Vector3D& gp) const final override;

  /// @copydoc SurfaceMaterial::material(size_t, size_t)
  ///
  /// @note the input parameter is ignored
  virtual const MaterialProperties*
  material(size_t ib0, size_t ib1) const final override;

  /// Output Method for std::ostream
  std::ostream&
  dump(std::ostream& sl) const final override;

private:
  /// The five different MaterialProperties
  MaterialProperties m_fullMaterial;
};

inline HomogeneousSurfaceMaterial*
HomogeneousSurfaceMaterial::clone() const
{
  return new HomogeneousSurfaceMaterial(*this);
}

inline const MaterialProperties*
HomogeneousSurfaceMaterial::material(const Vector2D&) const
{
  return (&m_fullMaterial);
}

inline const MaterialProperties*
HomogeneousSurfaceMaterial::material(const Vector3D&) const
{
  return (&m_fullMaterial);
}

inline const MaterialProperties*
    HomogeneousSurfaceMaterial::material(size_t, size_t) const
{
  return (&m_fullMaterial);
}

}  // end of namespace

#endif
