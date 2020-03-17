// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryStatics.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class LineBounds;

///  @class LineSurface
///
///  Base class for a linear surfaces in the TrackingGeometry
///  to describe dirft tube, straw like detectors or the Perigee
///  It inherits from Surface.
///
///  @note It leaves the type() method virtual, so it can not be instantiated
///
/// @image html LineSurface.png
class LineSurface : public Surface {
  friend Surface;

 protected:
  /// Constructor from Transform3D and bounds
  ///
  /// @param htrans The transform that positions the surface in the global frame
  /// @param radius The straw radius
  /// @param halez The half length in z
  LineSurface(std::shared_ptr<const Transform3D> htrans, double radius,
              double halez);

  /// Constructor from Transform3D and a shared bounds object
  ///
  /// @param htrans The transform that positions the surface in the global frame
  /// @param lbounds The bounds describing the straw dimensions, can be
  /// optionally nullptr
  LineSurface(std::shared_ptr<const Transform3D> htrans,
              std::shared_ptr<const LineBounds> lbounds = nullptr);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param lbounds The bounds describing the straw dimensions
  /// @param detelement for which this surface is (at least) one representation
  LineSurface(const std::shared_ptr<const LineBounds>& lbounds,
              const DetectorElementBase& detelement);

  /// Copy constructor
  ///
  /// @param other The source surface for copying
  LineSurface(const LineSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param transf is the additional transfrom applied after copying
  LineSurface(const GeometryContext& gctx, const LineSurface& other,
              const Transform3D& transf);

 public:
  /// Destructor - defaulted
  ~LineSurface() override = default;

  /// Default Constructor - deleted
  LineSurface() = delete;

  /// Assignment operator
  ///
  /// @param slsf is the source surface dor copying
  LineSurface& operator=(const LineSurface& other);

  /// Normal vector return
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position is ignored
  ///
  /// @return a Vector3D by value
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

  /// Return the measurement frame - this is needed for alignment, in particular
  ///
  /// for StraightLine and Perigee Surface
  ///  - the default implementation is the the RotationMatrix3D of the transform
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position where the measurement frame is
  /// constructed
  /// @param momentum is the momentum used for the measurement frame
  /// construction
  ///
  /// @return is a rotation matrix that indicates the measurement frame
  const RotationMatrix3D referenceFrame(const GeometryContext& gctx,
                                        const Vector3D& position,
                                        const Vector3D& momentum) const final;

  /// Initialize the jacobian from local to global
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param jacobian is the jacobian to be initialized
  /// @param position is the global position of the parameters
  /// @param direction is the direction at of the parameters
  ///
  /// @param pars is the paranmeters vector
  void initJacobianToGlobal(const GeometryContext& gctx,
                            BoundToFreeMatrix& jacobian,
                            const Vector3D& position, const Vector3D& direction,
                            const BoundVector& pars) const final;

  /// Calculate the form factors for the derivatives
  /// the calculation is identical for all surfaces where the
  /// reference frame does not depend on the direction
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the position of the paramters in global
  /// @param direction is the direction of the track
  /// @param rft is the transposed reference frame (avoids recalculation)
  /// @param jacobian is the transport jacobian
  ///
  /// @return a five-dim vector
  const BoundRowVector derivativeFactors(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction, const RotationMatrix3D& rft,
      const BoundToFreeMatrix& jacobian) const final;

  /// Local to global transformation
  /// for line surfaces the momentum is used in order to interpret the drift
  /// radius
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position to be transformed
  /// @param momentum is the global momentum (used to sign the closest approach)
  /// @param position is the global position which is filled
  void localToGlobal(const GeometryContext& gctx, const Vector2D& lposition,
                     const Vector3D& momentum, Vector3D& position) const final;

  /// Specified for LineSurface: global to local method without dynamic
  /// memory allocation
  ///
  /// This method is the true global->local transformation.<br>
  /// makes use of globalToLocal and indicates the sign of the Acts::eLOC_R by
  /// the given momentum
  ///
  /// The calculation of the sign of the radius (or \f$ d_0 \f$) can be done as
  /// follows:<br>
  /// May \f$ \vec d = \vec m - \vec c \f$ denote the difference between the
  /// center of the line and the global position of the measurement/predicted
  /// state,
  /// then \f$ \vec d
  /// \f$
  /// lies within the so
  /// called measurement plane.
  /// The measurement plane is determined by the two orthogonal vectors \f$
  /// \vec{measY}= \vec{Acts::eLOC_Z} \f$
  /// and \f$ \vec{measX} = \vec{measY} \times \frac{\vec{p}}{|\vec{p}|}
  /// \f$.<br>
  ///
  /// The sign of the radius (\f$ d_{0} \f$ ) is then defined by the projection
  /// of
  /// \f$ \vec{d} \f$
  /// onto \f$ \vec{measX} \f$:<br>
  /// \f$ sign = -sign(\vec{d} \cdot \vec{measX}) \f$
  ///
  /// \image html SignOfDriftCircleD0.gif
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
                     const Vector3D& momentum, Vector2D& lposition) const final;

  /// @brief Straight line intersection schema
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The global position as a starting point
  /// @param direction The global direction at the starting point
  ///        @note exptected to be normalized
  /// @param bcheck The boundary check directive for the estimate
  ///
  ///   <b>mathematical motivation:</b>
  ///   Given two lines in parameteric form:<br>
  ///   - @f$ \vec l_{a}(\lambda) = \vec m_a + \lambda \cdot \vec e_{a} @f$ <br>
  ///   - @f$ \vec l_{b}(\mu) = \vec m_b + \mu \cdot \vec e_{b} @f$ <br>
  ///   the vector between any two points on the two lines is given by:
  ///   - @f$ \vec s(\lambda, \mu) = \vec l_{b} - l_{a} = \vec m_{ab} + \mu
  ///   \cdot
  ///  \vec e_{b} - \lambda \cdot \vec e_{a} @f$, <br>
  ///   when @f$ \vec m_{ab} = \vec m_{b} - \vec m_{a} @f$.<br>
  ///   @f$ \vec s(u, \mu_0) @f$  denotes the vector between the two
  ///  closest points <br>
  ///   @f$ \vec l_{a,0} = l_{a}(u) @f$ and @f$ \vec l_{b,0} =
  ///  l_{b}(\mu_0) @f$ <br>
  ///   and is perpendicular to both, @f$ \vec e_{a} @f$ and @f$ \vec e_{b} @f$.
  ///
  ///   This results in a system of two linear equations:<br>
  ///   - (i) @f$ 0 = \vec s(u, \mu_0) \cdot \vec e_a = \vec m_ab \cdot
  ///  \vec e_a + \mu_0 \vec e_a \cdot \vec e_b - u @f$ <br>
  ///   - (ii) @f$ 0 = \vec s(u, \mu_0) \cdot \vec e_b = \vec m_ab \cdot
  ///  \vec e_b + \mu_0  - u \vec e_b \cdot \vec e_a @f$ <br>
  ///
  ///   Solving (i), (ii) for @f$ u @f$ and @f$ \mu_0 @f$ yields:
  ///   - @f$ u = \frac{(\vec m_ab \cdot \vec e_a)-(\vec m_ab \cdot \vec
  ///  e_b)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
  ///   - @f$ \mu_0 = - \frac{(\vec m_ab \cdot \vec e_b)-(\vec m_ab \cdot \vec
  ///  e_a)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
  ///
  /// @return is the intersection object
  Intersection intersectionEstimate(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction,
      const BoundaryCheck& bcheck = false) const final;

  /// the pathCorrection for derived classes with thickness
  /// is by definition 1 for LineSurfaces
  ///
  /// @note input parameters are ignored
  /// @note there's no material associated to the line surface
  double pathCorrection(const GeometryContext& gctx, const Vector3D& position,
                        const Vector3D& momentum) const override;

  /// This method returns the bounds of the Surface by reference */
  const SurfaceBounds& bounds() const final;

  /// Return properly formatted class name for screen output */
  std::string name() const override;

 protected:
  std::shared_ptr<const LineBounds> m_bounds;  ///< bounds (shared)

 private:
  /// helper function to apply the globalToLocal with out transform
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position
  /// @param momentum is the momentum
  /// @param lposition is the local position to be filled
  bool globalToLocalPlain(const GeometryContext& gctx, const Vector3D& position,
                          const Vector3D& momentum, Vector2D& lposition) const;
};

#include "Acts/Surfaces/detail/LineSurface.ipp"

}  // namespace Acts
