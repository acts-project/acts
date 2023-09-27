// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <ostream>
#include <string>

namespace Acts {

class DetectorElementBase;
class SurfaceBounds;
class ISurfaceMaterial;
class Layer;
class TrackingVolume;
class IVisualization3D;

/// Typedef of the surface intersection
using SurfaceIntersection = ObjectIntersection<Surface>;

/// @class Surface
///
/// @brief Abstract Base Class for tracking surfaces
///
/// The Surface class builds the core of the Acts Tracking Geometry.
/// All other geometrical objects are either extending the surface or
/// are built from it.
///
/// Surfaces are either owned by Detector elements or the Tracking Geometry,
/// in which case they are not copied within the data model objects.
///
class Surface : public virtual GeometryObject,
                public std::enable_shared_from_this<Surface> {
 public:
  /// @enum SurfaceType
  ///
  /// This enumerator simplifies the persistency & calculations,
  /// by saving a dynamic_cast, e.g. for persistency
  enum SurfaceType {
    Cone = 0,
    Cylinder = 1,
    Disc = 2,
    Perigee = 3,
    Plane = 4,
    Straw = 5,
    Curvilinear = 6,
    Other = 7
  };

 protected:
  /// Constructor with Transform3 as a shared object
  ///
  /// @param transform Transform3 positions the surface in 3D global space
  /// @note also acts as default constructor
  Surface(const Transform3& transform = Transform3::Identity());

  /// Copy constructor
  ///
  /// @note copy construction invalidates the association
  /// to detector element and layer
  ///
  /// @param other Source surface for copy.
  Surface(const Surface& other);

  /// Constructor fromt DetectorElementBase: Element proxy
  ///
  /// @param detelement Detector element which is represented by this surface
  Surface(const DetectorElementBase& detelement);

  /// Copy constructor with optional shift
  ///
  /// @note copy construction invalidates the association
  /// to detector element and layer
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other Source surface for copy
  /// @param shift Additional transform applied as: shift * transform
  Surface(const GeometryContext& gctx, const Surface& other,
          const Transform3& shift);

 public:
  virtual ~Surface();

  /// Factory for producing memory managed instances of Surface.
  /// Will forward all parameters and will attempt to find a suitable
  /// constructor.
  template <class T, typename... Args>
  static std::shared_ptr<T> makeShared(Args&&... args) {
    return std::shared_ptr<T>(new T(std::forward<Args>(args)...));
  }

  /// Retrieve a @c std::shared_ptr for this surface (non-const version)
  ///
  /// @note Will error if this was not created through the @c makeShared factory
  ///       since it needs access to the original reference. In C++14 this is
  ///       undefined behavior (but most likely implemented as a @c bad_weak_ptr
  ///       exception), in C++17 it is defined as that exception.
  /// @note Only call this if you need shared ownership of this object.
  ///
  /// @return The shared pointer
  std::shared_ptr<Surface> getSharedPtr();

  /// Retrieve a @c std::shared_ptr for this surface (const version)
  ///
  /// @note Will error if this was not created through the @c makeShared factory
  ///       since it needs access to the original reference. In C++14 this is
  ///       undefined behavior, but most likely implemented as a @c bad_weak_ptr
  ///       exception, in C++17 it is defined as that exception.
  /// @note Only call this if you need shared ownership of this object.
  ///
  /// @return The shared pointer
  std::shared_ptr<const Surface> getSharedPtr() const;

  /// Assignment operator
  /// @note copy construction invalidates the association
  /// to detector element and layer
  ///
  /// @param other Source surface for the assignment
  Surface& operator=(const Surface& other);

  /// Comparison (equality) operator
  /// The strategy for comparison is
  /// (a) first pointer comparison
  /// (b) then type comparison
  /// (c) then bounds comparison
  /// (d) then transform comparison
  ///
  /// @param other source surface for the comparison
  virtual bool operator==(const Surface& other) const;

  /// Comparison (non-equality) operator
  ///
  /// @param sf Source surface for the comparison
  virtual bool operator!=(const Surface& sf) const;

 public:
  /// Return method for the Surface type to avoid dynamic casts
  virtual SurfaceType type() const = 0;

  /// Return method for the surface Transform3 by reference
  /// In case a detector element is associated the surface transform
  /// is just forwarded to the detector element in order to keep the
  /// (mis-)alignment cache cetrally handled
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @return the contextual transform
  virtual const Transform3& transform(const GeometryContext& gctx) const;

  /// Return method for the surface center by reference
  /// @note the center is always recalculated in order to not keep a cache
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @return center position by value
  virtual Vector3 center(const GeometryContext& gctx) const;

  /// Return method for the normal vector of the surface
  /// The normal vector can only be generally defined at a given local position
  /// It requires a local position to be given (in general)
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position where the normal vector is
  /// constructed
  ///
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx,
                         const Vector2& lposition) const = 0;

  /// Return method for the normal vector of the surface
  /// The normal vector can only be generally defined at a given local position
  /// It requires a local position to be given (in general)
  ///
  /// @param position is the global position where the normal vector is
  /// constructed
  /// @param gctx The current geometry context object, e.g. alignment

  ///
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx,
                         const Vector3& position) const;

  /// Return method for the normal vector of the surface
  ///
  /// It will return a normal vector at the center() position
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  //
  /// @return normal vector by value
  virtual Vector3 normal(const GeometryContext& gctx) const {
    return normal(gctx, center(gctx));
  }

  /// Return method for SurfaceBounds
  /// @return SurfaceBounds by reference
  virtual const SurfaceBounds& bounds() const = 0;

  /// Return method for the associated Detector Element
  /// @return plain pointer to the DetectorElement, can be nullptr
  const DetectorElementBase* associatedDetectorElement() const;

  /// Return method for the associated Layer in which the surface is embedded
  /// @return Layer by plain pointer, can be nullptr
  const Layer* associatedLayer() const;

  /// Set Associated Layer
  /// Many surfaces can be associated to a Layer, but it might not be known yet
  /// during construction of the layer, this can be set afterwards
  ///
  /// @param lay the assignment Layer by reference
  void associateLayer(const Layer& lay);

  /// Return method for the associated Material to this surface
  /// @return SurfaceMaterial as plain pointer, can be nullptr
  const ISurfaceMaterial* surfaceMaterial() const;

  /// Return method for the shared pointer to the associated Material
  /// @return SurfaceMaterial as shared_pointer, can be nullptr
  const std::shared_ptr<const ISurfaceMaterial>& surfaceMaterialSharedPtr()
      const;

  /// Assign a detector element
  ///
  /// @param detelement Detector element which is represented by this surface
  void assignDetectorElement(const DetectorElementBase& detelement);

  /// Assign the surface material description
  ///
  /// The material is usually derived in a complicated way and loaded from
  /// a framework given source. As various surfaces may share the same source
  /// this is provided by a shared pointer
  ///
  /// @param material Material description associated to this surface
  void assignSurfaceMaterial(std::shared_ptr<const ISurfaceMaterial> material);

  /// The geometric onSurface method
  ///
  /// Geometrical check whether position is on Surface
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global position to be evaludated
  /// @param momentum global momentum (required for line-type surfaces)
  /// @param bcheck BoundaryCheck directive for this onSurface check
  ///
  /// @return boolean indication if operation was successful
  bool isOnSurface(const GeometryContext& gctx, const Vector3& position,
                   const Vector3& momentum,
                   const BoundaryCheck& bcheck = true) const;

  /// The insideBounds method for local positions
  ///
  /// @param lposition The local position to check
  /// @param bcheck BoundaryCheck directive for this onSurface check
  /// @return boolean indication if operation was successful
  virtual bool insideBounds(const Vector2& lposition,
                            const BoundaryCheck& bcheck = true) const;

  /// Local to global transformation
  /// Generalized local to global transformation for the surface types. Since
  /// some surface types need the global momentum/direction to resolve sign
  /// ambiguity this is also provided
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition local 2D position in specialized surface frame
  /// @param momentum global 3D momentum representation (optionally ignored)
  ///
  /// @return The global position by value
  virtual Vector3 localToGlobal(const GeometryContext& gctx,
                                const Vector2& lposition,
                                const Vector3& momentum) const = 0;

  /// Global to local transformation
  /// Generalized global to local transformation for the surface types. Since
  /// some surface types need the global momentum/direction to resolve sign
  /// ambiguity this is also provided
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param momentum global 3D momentum representation (optionally ignored)
  /// @param tolerance optional tolerance within which a point is considered
  /// valid on surface
  ///
  /// @return a Result<Vector2> which can be !ok() if the operation fails
  virtual Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& momentum,
      double tolerance = s_onSurfaceTolerance) const = 0;

  /// Return mehtod for the reference frame
  /// This is the frame in which the covariance matrix is defined (specialized
  /// by all surfaces)
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param momentum global 3D momentum representation (optionally ignored)
  ///
  /// @return RotationMatrix3 which defines the three axes of the measurement
  /// frame
  virtual Acts::RotationMatrix3 referenceFrame(const GeometryContext& gctx,
                                               const Vector3& position,
                                               const Vector3& momentum) const;

  /// Calculate the jacobian from local to global which the surface knows best,
  /// hence the calculation is done here.
  ///
  /// @note In principle, the input could also be a free parameters
  /// vector as it could be transformed to a bound parameters. But the transform
  /// might fail in case the parameters is not on surface. To avoid the check
  /// inside this function, it takes directly the bound parameters as input
  /// (then the check might be done where this function is called).
  ///
  /// @todo this mixes track parameterisation and geometry
  /// should move to :
  /// "Acts/EventData/detail/coordinate_transformations.hpp"
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param boundParams is the bound parameters vector
  ///
  /// @return Jacobian from local to global
  virtual BoundToFreeMatrix boundToFreeJacobian(
      const GeometryContext& gctx, const BoundVector& boundParams) const;

  /// Calculate the jacobian from global to local which the surface knows best,
  /// hence the calculation is done here.
  ///
  /// @note It assumes the input free parameters is on surface, hence no
  /// onSurface check is done inside this function.
  ///
  /// @todo this mixes track parameterisation and geometry
  /// should move to :
  /// "Acts/EventData/detail/coordinate_transformations.hpp"
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters is the free parameters
  ///
  /// @return Jacobian from global to local
  virtual FreeToBoundMatrix freeToBoundJacobian(
      const GeometryContext& gctx, const FreeVector& parameters) const;

  /// Calculate the derivative of path length at the geometry constraint or
  /// point-of-closest-approach w.r.t. free parameters. The calculation is
  /// identical for all surfaces where the reference frame does not depend on
  /// the direction
  ///
  /// @todo this mixes track parameterisation and geometry
  /// should move to :
  /// "Acts/EventData/detail/coordinate_transformations.hpp"
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters is the free parameters
  ///
  /// @return Derivative of path length w.r.t. free parameters
  virtual FreeToPathMatrix freeToPathDerivative(
      const GeometryContext& gctx, const FreeVector& parameters) const;

  /// Calucation of the path correction for incident
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param direction global 3D momentum direction
  ///
  /// @return Path correction with respect to the nominal incident.
  virtual double pathCorrection(const GeometryContext& gctx,
                                const Vector3& position,
                                const Vector3& direction) const = 0;

  /// Straight line intersection schema from position/direction
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to start from
  /// @param direction The direction at start
  /// @param bcheck the Boundary Check
  ///
  /// @return SurfaceIntersection object (contains intersection & surface)
  virtual SurfaceIntersection intersect(const GeometryContext& gctx,
                                        const Vector3& position,
                                        const Vector3& direction,
                                        const BoundaryCheck& bcheck) const = 0;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param sl is the ostream to be dumped into
  virtual std::ostream& toStream(const GeometryContext& gctx,
                                 std::ostream& sl) const;

  /// Output into a std::string
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  std::string toString(const GeometryContext& gctx) const;

  /// Return properly formatted class name
  virtual std::string name() const = 0;

  /// Return a Polyhedron for this object
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, if the lseg
  /// is set to one, only the corners and the extrema are given,
  /// otherwise it represents the number of segments for a full 2*M_PI
  /// circle and is scaled to the relevant sector
  ///
  /// @note An internal surface transform can invalidate the extrema
  /// in the transformed space
  ///
  /// @return A list of vertices and a face/facett description of it
  virtual Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                              size_t lseg) const = 0;

  /// The derivative of bound track parameters w.r.t. alignment
  /// parameters of its reference surface (i.e. local frame origin in
  /// global 3D Cartesian coordinates and its rotation represented with
  /// extrinsic Euler angles)
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// change of alignment parameters
  /// @param parameters is the free parameters
  /// @param pathDerivative is the derivative of free parameters w.r.t. path
  /// length
  ///
  /// @return Derivative of bound track parameters w.r.t. local frame
  /// alignment parameters
  AlignmentToBoundMatrix alignmentToBoundDerivative(
      const GeometryContext& gctx, const FreeVector& parameters,
      const FreeVector& pathDerivative) const;

  /// Calculate the derivative of path length at the geometry constraint or
  /// point-of-closest-approach w.r.t. alignment parameters of the surface (i.e.
  /// local frame origin in global 3D Cartesian coordinates and its rotation
  /// represented with extrinsic Euler angles)
  ///
  /// @note Re-implementation is needed for surface whose intersection with
  /// track is not its local xy plane, e.g. LineSurface, CylinderSurface and
  /// ConeSurface
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters is the free parameters
  ///
  /// @return Derivative of path length w.r.t. the alignment parameters
  virtual AlignmentToPathMatrix alignmentToPathDerivative(
      const GeometryContext& gctx, const FreeVector& parameters) const;

  /// Calculate the derivative of bound track parameters local position w.r.t.
  /// position in local 3D Cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position of the paramters in global
  ///
  /// @return Derivative of bound local position w.r.t. position in local 3D
  /// cartesian coordinates
  virtual ActsMatrix<2, 3> localCartesianToBoundLocalDerivative(
      const GeometryContext& gctx, const Vector3& position) const = 0;

 protected:
  /// Transform3 definition that positions
  /// (translation, rotation) the surface in global space
  Transform3 m_transform = Transform3::Identity();

  /// Pointer to the a DetectorElementBase
  const DetectorElementBase* m_associatedDetElement{nullptr};

  /// The associated layer Layer - layer in which the Surface is be embedded,
  /// nullptr if not associated
  const Layer* m_associatedLayer{nullptr};

  /// The assoicated TrackingVolume - tracking volume in case the surface is a
  /// boundary surface, nullptr if not associated
  const TrackingVolume* m_associatedTrackingVolume{nullptr};

  /// Possibility to attach a material descrption
  std::shared_ptr<const ISurfaceMaterial> m_surfaceMaterial;

 private:
  /// Calculate the derivative of bound track parameters w.r.t.
  /// alignment parameters of its reference surface (i.e. origin in global 3D
  /// Cartesian coordinates and its rotation represented with extrinsic Euler
  /// angles) without any path correction
  ///
  /// @note This function should be used together with alignment to path
  /// derivative to get the full alignment to bound derivatives
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters is the free parameters
  ///
  /// @return Derivative of bound track parameters w.r.t. local frame alignment
  /// parameters without path correction
  AlignmentToBoundMatrix alignmentToBoundDerivativeWithoutCorrection(
      const GeometryContext& gctx, const FreeVector& parameters) const;
};

/// Print surface information to the provided stream. Internally invokes the
/// `surface.toStream(...)`-method. This can be easily used e.g. like `std::cout
/// << std::tie(surface, geometryContext);`
inline std::ostream& operator<<(
    std::ostream& os,
    const std::tuple<const Surface&, const GeometryContext&>& tup) {
  const auto [surface, gctx] = tup;
  surface.toStream(gctx, os);
  return os;
}

}  // namespace Acts
