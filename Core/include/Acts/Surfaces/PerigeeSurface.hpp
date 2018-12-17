// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// PerigeeSurface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryStatics.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

/// @class PerigeeSurface
///
/// Class describing the Line to which the Perigee refers to.
/// The Surface axis is fixed to be the z-axis of the Tracking frame.
/// It inherits from StraingLineSurface.
///
/// @image html LineSurface.png
class PerigeeSurface : public LineSurface
{
  friend Surface;

  PerigeeSurface() = delete;

protected:
  /// Constructor from GlobalPosition
  ///
  /// @param gp position where the perigee is centered
  PerigeeSurface(const Vector3D& gp);

  /// Constructor with a Transform - needed for tilt
  ///
  /// @param tTransform is the transform for position and tilting
  PerigeeSurface(std::shared_ptr<const Transform3D> tTransform);

  /// Copy constructor
  ///
  /// @param other is the source surface to be copied
  PerigeeSurface(const PerigeeSurface& other);

  /// Copy constructor with shift
  ///
  /// @param other is the source surface to be copied
  /// @param shift is the transformed applied after copying
  PerigeeSurface(const PerigeeSurface& other, const Transform3D& shift);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  PerigeeSurface(const variant_data& vardata);

public:
  /// Destructor - defaulted
  ~PerigeeSurface() override = default;

  /// Clone method. Uses the copy constructor a new position can optionally be
  /// given a shift.
  ///
  /// @param shift additional, optional shift
  std::shared_ptr<PerigeeSurface>
  clone(const Transform3D* shift = nullptr) const;

  /// Assignment operator
  ///
  /// @param other is the source surface to be assigned
  PerigeeSurface&
  operator=(const PerigeeSurface& other);

  /// Return the surface type
  SurfaceType
  type() const final;

  /// Return properly formatted class name for screen output */
  std::string
  name() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream&
  dump(std::ostream& sl) const final;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

private:
  /// Clone method. Uses the copy constructor a new position can optionally be
  /// given a shift.
  ///
  /// @param shift additional, optional shift
  PerigeeSurface*
  clone_impl(const Transform3D* shift = nullptr) const override;
};

}  // namespace