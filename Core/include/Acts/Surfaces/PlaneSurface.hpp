// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryStatics.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/Definitions.hpp"

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
class PlaneSurface : public Surface {
  friend Surface;

 protected:
  /// Default Constructor - needed for persistency
  PlaneSurface();

  /// Copy Constructor
  ///
  /// @param psf is the source surface for the copy
  PlaneSurface(const PlaneSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param transf is the additional transfrom applied after copying
  PlaneSurface(const GeometryContext& gctx, const PlaneSurface& other,
               const Transform3D& transf);

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
               const DetectorElementBase& detelement);

  /// Constructor for Planes with (optional) shared bounds object
  ///
  /// @param htrans transform in 3D that positions this surface
  /// @param pbounds bounds object to describe the actual surface area
  PlaneSurface(std::shared_ptr<const Transform3D> htrans,
               std::shared_ptr<const PlanarBounds> pbounds = nullptr);

 public:
  /// Destructor - defaulted
  ~PlaneSurface() override = default;

  /// Assignment operator
  ///
  /// @param other The source PlaneSurface for assignment
  PlaneSurface& operator=(const PlaneSurface& other);

  /// Clone method into a concrete type of PlaneSurface with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  std::shared_ptr<PlaneSurface> clone(const GeometryContext& gctx,
                                      const Transform3D& shift) const;

  /// Normal vector return
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position is ignored
  ///
  /// return a Vector3D by value
  const Vector3D normal(const GeometryContext& gctx,
                        const Vector2D& lposition) const final;

  /// Normal vector return without argument
  using Surface::normal;

  /// The binning position is the position calcualted
  /// for a certain binning type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the binning type to be used
  ///
  /// @return position that can beused for this binning
  const Vector3D binningPosition(const GeometryContext& gctx,
                                 BinningValue bValue) const final;

  /// Return the surface type
  SurfaceType type() const override;

  /// Return method for bounds object of this surfrace
  const SurfaceBounds& bounds() const override;

  /// Local to global transformation
  /// For planar surfaces the momentum is ignroed in the local to global
  /// transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition local 2D position in specialized surface frame
  /// @param momentum global 3D momentum representation (optionally ignored)
  /// @param position global 3D position to be filled (given by reference for
  /// method symmetry)
  void localToGlobal(const GeometryContext& gctx, const Vector2D& lposition,
                     const Vector3D& momentum,
                     Vector3D& position) const override;

  /// Global to local transformation
  /// For planar surfaces the momentum is ignroed in the global to local
  /// transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param momentum global 3D momentum representation (optionally ignored)
  /// @param lposition local 2D position to be filled (given by reference for
  /// method symmetry)
  ///
  /// @return boolean indication if operation was successful (fail means global
  /// position was not on surface)
  bool globalToLocal(const GeometryContext& gctx, const Vector3D& position,
                     const Vector3D& momentum,
                     Vector2D& lposition) const override;

  /// Method that calculates the correction due to incident angle
  ///
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param direction global 3D momentum direction (ignored for PlaneSurface)
  /// @note this is the final implementation of the pathCorrection function
  ///
  /// @return a double representing the scaling factor
  double pathCorrection(const GeometryContext& gctx, const Vector3D& position,
                        const Vector3D& direction) const final;

  /// @brief Straight line intersection schema
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The start position of the intersection attempt
  /// @param direction The direction of the interesection attempt,
  /// (@note expected to be normalized)
  /// @param bcheck The boundary check directive
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
  Intersection intersectionEstimate(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction,
      const BoundaryCheck& bcheck = false) const final;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, it represents
  /// the full 2*M_PI coverange, if lseg is set to 1 only the extrema
  /// are given
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      size_t lseg) const override;

  /// Return properly formatted class name for screen output
  std::string name() const override;

 protected:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;

 private:
  /// Clone method implementation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  PlaneSurface* clone_impl(const GeometryContext& gctx,
                           const Transform3D& shift) const override;
};

#include "Acts/Surfaces/detail/PlaneSurface.ipp"

}  // end of namespace Acts
