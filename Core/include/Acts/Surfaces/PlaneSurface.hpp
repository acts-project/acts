// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <limits>
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryStatics.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

class DetectorElementBase;

/// @class PlaneSurface
///
/// Class for a planaer in the TrackingGeometry.
///
/// The PlaneSurface extends the Surface class with the possibility to
/// convert local to global positions (vice versa).
///
/// @image html PlaneSurface.png
///
class PlaneSurface : public Surface
{
  friend Surface;

protected:
  /// Default Constructor - needed for persistency
  PlaneSurface();

  /// Copy Constructor
  ///
  /// @param psf is the source surface for the copy
  PlaneSurface(const PlaneSurface& other);

  /// Copy Constructor with shift
  ///
  /// @param other The source surface for the copy
  /// @param transf The transformation that positions the surface in space
  PlaneSurface(const PlaneSurface& other, const Transform3D& transf);

  /// Dedicated Constructor with normal vector
  /// This is for curvilinear surfaces which are by definition boundless
  ///
  /// @param center is the center position of the surface
  /// @param normal is thenormal vector of the plane surface
  PlaneSurface(const Vector3D& center, const Vector3D& normal);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param pbounds are the provided planar bounds (shared)
  /// @param detelement is the linked detector element to this surface
  PlaneSurface(const std::shared_ptr<const PlanarBounds>& pbounds,
               const DetectorElementBase&                 detelement);

  /// Constructor for Planes with (optional) shared bounds object
  ///
  /// @param htrans transform in 3D that positions this surface
  /// @param pbounds bounds object to describe the actual surface area
  PlaneSurface(std::shared_ptr<const Transform3D>  htrans,
               std::shared_ptr<const PlanarBounds> pbounds = nullptr);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  PlaneSurface(const variant_data& vardata);

public:
  /// Destructor - defaulted
  ~PlaneSurface() override = default;

  /// Assignment operator
  ///
  /// @param other The source PlaneSurface for assignment
  PlaneSurface&
  operator=(const PlaneSurface& other);

  /// Clone method. Uses the copy constructor a new position can optionally be
  /// given a shift.
  ///
  /// @param shift additional, optional shift
  std::shared_ptr<PlaneSurface>
  clone(const Transform3D* shift = nullptr) const;

  /// Normal vector return
  ///
  /// @param lpos is the local position is ignored
  /// return a Vector3D by value
  const Vector3D
  normal(const Vector2D& lpos) const final;

  /// Normal vector return without argument
  using Surface::normal;

  /// The binning position is the position calcualted
  /// for a certain binning type
  ///
  /// @param bValue is the binning type to be used
  ///
  /// @return position that can beused for this binning
  const Vector3D
  binningPosition(BinningValue bValue) const final;

  /// Return the surface type
  SurfaceType
  type() const override;

  /// Return method for bounds object of this surfrace
  const SurfaceBounds&
  bounds() const override;

  /// Local to global transformation
  /// For planar surfaces the momentum is ignroed in the local to global
  /// transformation
  ///
  /// @param lpos local 2D posittion in specialized surface frame
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param gpos global 3D position to be filled (given by reference for method
  /// symmetry)
  void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const override;

  /// Global to local transformation
  /// For planar surfaces the momentum is ignroed in the global to local
  /// transformation
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param lpos local 2D position to be filled (given by reference for method
  /// symmetry)
  ///
  /// @return boolean indication if operation was successful (fail means global
  /// position was not on surface)
  bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const override;

  /// Method that calculates the correction due to incident angle
  ///
  /// @param pos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param mom global 3D momentum representation (optionally ignored)
  ///
  /// @note this is the final implementation of the pathCorrection function
  ///
  /// @return a double representing the scaling factor
  double
  pathCorrection(const Vector3D& pos, const Vector3D& mom) const final;

  /// @brief Fast straight line intersection schema
  ///
  /// @param gpos The start position of the intersection attempt
  /// @param gdir The direction of the interesection attempt,
  ///       @note expected to be normalized
  /// @param navDir The navigation direction with respect to the momentum
  /// @param bcheck The boundary check directive
  /// @param correct is a corrector function (e.g. for curvature correction)
  ///
  /// <b>mathematical motivation:</b>
  ///
  /// the equation of the plane is given by: <br>
  /// @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  /// where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane,  @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point
  /// on the plane and @f$ \vec x = (x,y,z) @f$ all possible points
  /// on the plane.<br>
  ///
  /// Given a line with:<br>
  /// @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  /// the solution for @f$ u @f$ can be written:
  /// @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  /// If the denominator is 0 then the line lies:
  /// - either in the plane
  /// - perpendicular to the normal of the plane
  ///
  /// @return the Intersection object
  Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gdir,
                       NavigationDirection  navDir  = forward,
                       const BoundaryCheck& bcheck  = false,
                       CorrFnc              correct = nullptr) const final;

  /// Return properly formatted class name for screen output
  std::string
  name() const override;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

protected:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;

private:
  /// Clone method. Uses the copy constructor a new position can optionally be
  /// given a shift.
  ///
  /// @param shift additional, optional shift
  PlaneSurface*
  clone_impl(const Transform3D* shift = nullptr) const override;
};

#include "Acts/Surfaces/detail/PlaneSurface.ipp"

}  // end of namespace Acts