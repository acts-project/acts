// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscSurface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Identifier.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

class DetectorElementBase;

/// @class DiscSurface
///
/// Class for a DiscSurface in the, it inherits from Surface.
///
/// The DiscSurface has a polar local coordinate system, with
/// (r,phi) describing the coordinates.
///
/// The surface transform positions the disc such, that the origin
/// is at r=0, independent of the provided DiscBounds. The z-axis
/// is the normal vector of the Disc, being perpendicular to the
/// radial direction.
///
/// The disc surface is the only surface type for which the
/// cavairance matrix is NOT given in the reference frame.
/// A conversion from polar to cartesian coordinates needs
/// to happen to transfer the local coordinates onto the
/// cartesian reference frame coordinates.
///
/// @image html DiscSurface.png
///
class DiscSurface : public Surface
{
public:
  /// Default Constructor is deleted
  DiscSurface() = delete;

  /// Constructor for Discs from Transform3D, \f$ r_{min}, r_{max} \f$
  ///
  /// @param htrans is transform that places the disc in the global 3D space
  /// (can be nullptr)
  /// @param rmin is the inner radius of the disc surface
  /// @param rmax is the outer radius of the disc surface
  /// @param hphisec is the opening angle of the disc surface and is optional
  ///        the default is a full disc
  DiscSurface(std::shared_ptr<const Transform3D> htrans,
              double                             rmin,
              double                             rmax,
              double                             hphisec = M_PI);

  /// Constructor for Discs from Transform3D, \f$ r_{min}, r_{max}, hx_{min},
  /// hx_{max} \f$
  /// This is n this case you have DiscTrapezoidalBounds
  ///
  /// @param htrans is transform that places the disc in the global 3D space
  /// (can be nullptr)
  /// @param minhalfx is the half length in x at minimal r
  /// @param maxhalfx is the half length in x at maximal r
  /// @param rmin is the inner radius of the disc surface
  /// @param rmax is the outer radius of the disc surface
  /// @param avephi is the position in phi (default is 0.)
  /// @param stereo is the optional stereo angle
  DiscSurface(std::shared_ptr<const Transform3D> htrans,
              double                             minhalfx,
              double                             maxhalfx,
              double                             rmin,
              double                             rmax,
              double                             avephi = 0.,
              double                             stereo = 0.);

  /// Constructor for Discs from Transform3D and shared DiscBounds
  ///
  /// @param htrans is the transform that positions the disc in the global 3D
  /// frame
  /// @param dbounds are the disc bounds describing the surface coverage
  DiscSurface(std::shared_ptr<const Transform3D> htrans,
              std::shared_ptr<const DiscBounds>  dbounds = nullptr);

  /// Constructor from detector element and identifier
  /// @note the surface only acts as a proxy of the detector element
  ///
  /// @param dbounds are the disc bounds associated to this surface, must not be
  /// nullptr
  /// @param detelement is the detector element that is represented by this
  /// surface
  /// @param identifier is the optional identifier in case one detector element
  /// owns more than 1 surface
  DiscSurface(std::shared_ptr<const DiscBounds> dbounds,
              const DetectorElementBase&        detelement,
              const Identifier&                 identifier = Identifier());

  /// Copy Constructor
  ///
  /// @param other is the source surface for the copy
  DiscSurface(const DiscSurface& other);

  /// Copy Constructor with shift
  ///
  /// @param other is the source sourface for the copy
  /// @param transf is the additional transform applied to the surface
  DiscSurface(const DiscSurface& other, const Transform3D& transf);

  /// Constructor which accepts @c variant_data
  ///
  /// @param data the @c variant_data to build from
  DiscSurface(const variant_data& data);

  virtual ~DiscSurface();

  /// Assignement operator
  ///
  /// @param other is the source sourface for the assignment
  DiscSurface&
  operator=(const DiscSurface& other);

  /// Virtual constructor - shift can be given optionally
  ///
  /// @param shift the otional transform applied after cloning
  virtual DiscSurface*
  clone(const Transform3D* shift = nullptr) const final override;

  /// Return the surface type
  virtual SurfaceType
  type() const override;

  /// Normal vector return
  ///
  /// @param lpos is the local position is ignored
  /// return a Vector3D by value
  const Vector3D
  normal(const Vector2D& lpos) const override final;

  /// Normal vector return
  ///
  /// @param lpos is the global position is ignored
  /// return a Vector3D by value
  const Vector3D
  normal(const Vector3D& gpos) const override final;

  /// Normal vector return
  ///
  /// @note No param, this overload resolves default parameter ambiguity
  const Vector3D
  normal() const;

  /// The binning position is the position calcualted
  /// for a certain binning type
  ///
  /// @param bValue is the binning type to be used
  /// @return position that can beused for this binning
  virtual const Vector3D
  binningPosition(BinningValue bValue) const final;

  /// This method returns the bounds by reference
  const SurfaceBounds&
  bounds() const final override;

  /// This method returns true if the GlobalPosition is on the Surface for both,
  /// within or without check of whether the local position is inside boundaries
  /// or not
  ///
  /// @param gpos is the global position to be checked
  /// @param bcheck is the boundary check directive
  /// @return bollean that indicates if the position is on surface
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bcheck = true) const final override;

  /// Local to global transformation
  /// For planar surfaces the momentum is ignroed in the local to global
  /// transformation
  ///
  /// @param lpos local 2D posittion in specialized surface frame
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param gpos global 3D position to be filled (given by reference for method
  /// symmetry)
  ///
  /// @note the momentum is ignored for Disc surfaces in this calculateion
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const final override;

  /// Global to local transformation
  /// @note the momentum is ignored for Disc surfaces in this calculateion
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param lpos local 2D position to be filled (given by reference for method
  /// symmetry)
  /// @return boolean indication if operation was successful (fail means global
  /// position was not on surface)
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const final override;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lpolar is a local position in polar coordinates
  /// @return values is local 2D position in cartesian coordinates  @todo check
  const Vector2D
  localPolarToCartesian(const Vector2D& lpolar) const;

  /// Special method for Disc surface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lcart is local 2D position in cartesian coordinates
  /// @return value is a local position in polar coordinates
  const Vector2D
  localCartesianToPolar(const Vector2D& lcart) const;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lpolar is a local position in polar coordinates
  /// @return values is local 2D position in cartesian coordinates
  const Vector2D
  localPolarToLocalCartesian(const Vector2D& lpolar) const;

  /// Special method for DiscSurface :  local<->global transformation when
  /// provided cartesian coordinates
  ///
  /// @param lcart is local 2D position in cartesian coordinates
  /// @return value is a global cartesian 3D position
  const Vector3D
  localCartesianToGlobal(const Vector2D& lcart) const;

  /// Special method for DiscSurface : global<->local from cartesian coordinates
  ///
  /// @param gpos is a global cartesian 3D position
  /// @param tol is the absoltue tolerance parameter
  /// @return value is a local polar
  const Vector2D
  globalToLocalCartesian(const Vector3D& gpos, double tol = 0.) const;

  /// Initialize the jacobian from local to global
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param jac is the jacobian to be initialized
  /// @param pos is the global position of the parameters
  /// @param dir is the direction at of the parameters
  /// @param pars is the paranmeters vector
  virtual void
      initJacobianToGlobal(ActsMatrixD<7, 5>& jac,
                           const Vector3D&       gpos,
                           const Vector3D&       dir,
                           const ActsVectorD<5>& pars) const final override;

  /// Initialize the jacobian from global to local
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param jac is the jacobian to be initialized
  /// @param pos is the global position of the parameters
  /// @param dir is the direction at of the parameters
  /// @param pars is the parameter vector
  ///
  /// @return the transposed reference frame (avoids recalculation)
  virtual const RotationMatrix3D
      initJacobianToLocal(ActsMatrixD<5, 7>& jac,
                          const Vector3D& gpos,
                          const Vector3D& dir) const final override;

  /// Path correction due to incident of the track
  ///
  /// @param gpos is the global position as a starting point
  /// @param mom is the global momentum at the starting point
  /// @return is the correction factor due to incident
  double
  pathCorrection(const Vector3D& gpos,
                 const Vector3D& mom) const final override;

  /// @brief Fast straight line intersection schema
  ///
  /// navDir=anyDirection is to provide the closest solution
  ///
  /// @param gpos is the global position as a starting point
  /// @param gdir is the global direction at the starting point
  ///        @note expected to be normalized (no checking)
  /// @param navDir is a navigation direction
  /// @param bcheck is the boundary check prescription
  /// @param correct is a corrector function (e.g. for curvature correction)
  ///
  ///  <b>mathematical motivation:</b>
  ///
  /// the equation of the plane is given by: <br>
  /// @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  /// where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane, @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on
  /// the plane and @f$ \vec x = (x,y,z) @f$ all possible points
  /// on the plane.<br>
  /// Given a line with:<br>
  /// @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  /// the solution for @f$ u @f$ can be written:
  /// @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  /// If the denominator is 0 then the line lies:
  /// - either in the plane
  /// - perpendicular to the normal of the plane
  ///
  /// @return is the surface intersection object
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gdir,
                       NavigationDirection  navDir = forward,
                       const BoundaryCheck& bcheck = false,
                       CorrFnc correct = nullptr) const final override;

  /// Return properly formatted class name for screen output
  virtual std::string
  name() const override;

  /// Return a PolyhedronRepresentation for this object
  /// @param l0div Number of divisions along l0 (phi)
  /// @param l1div Number of divisions along l1 (r)
  virtual PolyhedronRepresentation
  polyhedronRepresentation(size_t l0div = 10, size_t l1div = 1) const;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  virtual variant_data
  toVariantData() const override;

protected:
  std::shared_ptr<const DiscBounds> m_bounds;  ///< bounds (shared)
};

#include "detail/DiscSurface.ipp"

}  // end of namespace
