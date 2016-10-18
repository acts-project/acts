// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// PerigeeSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_PERIGEESURFACE_H
#define ACTS_SURFACES_PERIGEESURFACE_H 1

#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Surfaces/LineSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometryStatics.hpp"

namespace Acts {

/// @class PerigeeSurface
///
/// Class describing the Line to which the Perigee refers to.
/// The Surface axis is fixed to be the z-axis of the Tracking frame.
/// It inherits from StraingLineSurface.

class PerigeeSurface : public LineSurface
{
public:
  /// Default Constructor - deleted */
  PerigeeSurface() = delete;

  /// Constructor from GlobalPosition
  ///
  /// @param gpos position where the perigee is centered
  PerigeeSurface(const Vector3D& gpos);

  /// Constructor with a Transform - needed for tilt
  ///
  /// @param tTransform is the transform for position and tilting
  PerigeeSurface(std::shared_ptr<Transform3D> tTransform);

  /// Copy constructor
  ///
  /// @param pesf is the source surface to be copied
  PerigeeSurface(const PerigeeSurface& pesf);

  /// Copy constructor with shift
  ///
  /// @param pesf is the source surface to be copied
  /// @param transf is the transformed applied after copying
  PerigeeSurface(const PerigeeSurface& pesf, const Transform3D& transf);

  /// Destructor
  virtual ~PerigeeSurface();

  /// Virtual constructor
  ///
  /// @param shift is the potential shift that is applied after cloning
  virtual PerigeeSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /// Assignment operator
  ///
  /// @param pesf is the source surface to be assigned
  PerigeeSurface&
  operator=(const PerigeeSurface& pesf);

  /// Return the surface type
  virtual SurfaceType
  type() const override
  {
    return Surface::Perigee;
  }

  /// Return properly formatted class name for screen output */
  virtual std::string
  name() const override
  {
    return "Acts::PerigeeSurface";
  }

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  virtual std::ostream&
  dump(std::ostream& sl) const override;
};

inline PerigeeSurface*
PerigeeSurface::clone(const Transform3D* shift) const
{
  if (shift) return new PerigeeSurface(*this, *shift);
  return new PerigeeSurface(*this);
}

}  // end of namespace

#endif  // ACTS_SURFACESPERIGEESURFACE_H
