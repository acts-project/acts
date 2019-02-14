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

  /// Copy constructor - with shift
  ///
  /// @param ctx Is the payload/context object to be used for
  ///        delegating the event or thread context
  /// @param other is the source cone surface
  /// @param transf is the additional transfrom applied after copying
  PerigeeSurface(Context               ctx,
                 const PerigeeSurface& other,
                 const Transform3D&    transf);

public:
  /// Destructor - defaulted
  ~PerigeeSurface() override = default;

  /// Clone method into a concrete type of PerigeeSurface with shift
  ///
  /// @param ctx Is the payload/context object to be used for
  ///        delegating the event or thread context
  /// @param shift applied to the surface
  std::shared_ptr<PerigeeSurface>
  clone(Context ctx, const Transform3D& shift) const;

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
  /// @param ctx Is the payload/context object to be used for
  ///        delegating the event or thread context
  /// @param sl is the ostream to be dumped into
  ///
  /// @return ostreamn obect which was streamed into
  std::ostream&
  toStream(Context ctx, std::ostream& sl) const final;

private:
  /// Clone method implementation
  ///
  /// @param ctx Is the payload/context object to be used for
  ///        delegating the event or thread context
  /// @param shift applied to the surface
  PerigeeSurface*
  clone_impl(Context ctx, const Transform3D& shift) const override;
};

}  // namespace Acts
