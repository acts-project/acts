// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Surface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Detector/DetectorElementBase.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryObject.hpp"
#include "Acts/Utilities/GeometryStatics.hpp"
#include "Acts/Utilities/Identifier.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

class DetectorElementBase;
class SurfaceBounds;
class SurfaceMaterial;
class Layer;
class TrackingVolume;

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
class Surface : public virtual GeometryObject
{
public:
  /// @enum SurfaceType
  ///
  /// This enumerator simplifies the persistency & calculations,
  /// by saving a dynamic_cast, e.g. for persistency
  ///
  enum SurfaceType {
    Cone        = 0,
    Cylinder    = 1,
    Disc        = 2,
    Perigee     = 3,
    Plane       = 4,
    Straw       = 5,
    Curvilinear = 6,
    Other       = 7
  };

  /// Typedef of the surface intersection
  typedef ObjectIntersection<Surface> SurfaceIntersection;

  /// Constructor with Transform3D as a shared object
  ///
  /// @param htrans Transform3D positions the surface in 3D global space
  /// @note also acts as default constructor
  Surface(std::shared_ptr<const Transform3D> htrans = nullptr);

  /// Copy constructor
  /// @note copy construction invalidates the association
  /// to detector element and identifier
  ///
  /// @param other Source surface for copy.
  Surface(const Surface& other);

  /// Copy constructor with shift
  /// @note copy construction invalidates the association
  /// to detector element and identifier
  ///
  /// @param other Source surface for copy
  /// @param transf Additional transform applied after copying from the source
  Surface(const Surface& other, const Transform3D& transf);

  /// Constructor fromt DetectorElementBase and (optional) Identifier
  ///
  /// @param detelement Detector element which is represented by this surface
  /// @param id Optional identifier if more than one surface are associated to a
  /// detector element
  Surface(const DetectorElementBase& detelement,
          const Identifier&          id = Identifier());

  virtual ~Surface();

  /// Assignment operator is not allowed
  /// @note copy construction invalidates the association
  /// to detector element and identifier
  ///
  /// @param other Source surface for the assignment
  Surface&
  operator=(const Surface& other);

  /// Comparison (equality) operator
  /// The strategy for comparison is
  /// (a) first pointer comparison
  /// (b) then type comparison
  /// (c) then bounds comparison
  /// (d) then transform comparison
  ///
  /// @param other source surface for the comparison
  virtual bool
  operator==(const Surface& other) const;

  /// Comparison (non-equality) operator
  ///
  /// @param other Source surface for the comparison
  virtual bool
  operator!=(const Surface& other) const;

  /// Implicit constructor
  /// uses the copy constructor a new position can optionally be given
  ///
  /// @param shift is an additional (optional shift can be given)
  /// shift after cloning
  virtual Surface*
  clone(const Transform3D* shift = nullptr) const = 0;

  /// Return method for the Surface type to avoid dynamic casts
  virtual SurfaceType
  type() const = 0;

  /// Return method for the surface Transform3D by reference
  /// In case a detector element is associated the surface transform
  /// is just forwarded to the detector element in order to keep the
  /// (mis-)alignment cache cetrally handled
  virtual const Transform3D&
  transform() const;

  /// Return method for the surface center by reference
  /// @note the center is always recalculated in order to not keep a cache
  /// @return center position by value
  virtual const Vector3D
  center() const;

  /// Return method for the normal vector of the surface
  /// The normal vector can only be generally defined at a given local position
  /// It requires a local position to be given (in general)
  ///
  /// @param lpos is the local position where the normal verctor is constructed
  /// @return normal vector by value
  virtual const Vector3D
  normal(const Vector2D& lpos) const = 0;

  /// Return method for the normal vector of the surface
  /// The normal vector can only be generally defined at a given local position
  /// It requires a local position to be given (in general)
  ///
  /// @param pos is the global position where the normal vector is constructed
  /// @return normal vector by value
  virtual const Vector3D
  normal(const Vector3D& pos) const;

  /// Return method for SurfaceBounds
  /// @return SurfaceBounds by reference
  virtual const SurfaceBounds&
  bounds() const = 0;

  /// Return method for the associated Detector Element
  /// @return plain pointer to the DetectorElement, can be nullptr
  const DetectorElementBase*
  associatedDetectorElement() const;

  /// Return method for the associated Identifier
  /// @return Identifier by value, can be invalid
  const Identifier
  associatedIdentifier() const;

  /// Return method for the associated Layer in which the surface is embedded
  /// @return Layer by plain pointer, can be nullptr
  const Layer*
  associatedLayer() const;

  /// Return method for the associated Material to this surface
  /// @return SurfaceMaterial as plain pointer, can be nullptr
  const SurfaceMaterial*
  associatedMaterial() const;

  /// Set Associated Layer
  /// Many surfaces can be associated to a Layer, but it might not be known yet
  /// during construction of the layer, this can be set afterwards
  ///
  /// @param lay the assignment Layer by reference
  void
  associateLayer(const Layer& lay);

  /// Set Associated SurfaceMaterial
  /// The material is usually derived in a complicated way and loaded from
  /// a framework given source. As various srufaces may share the same
  ///
  /// @param material Material description this given and stored as a shared
  /// pointer
  void
  setAssociatedMaterial(std::shared_ptr<const SurfaceMaterial> material);

  /// The templated Parameters onSurface method
  /// In order to avoid unneccessary geometrical operations, it checks on the
  /// surface pointer first. If that check fails, it calls the geometrical
  /// check isOnSurface
  ///
  /// @tparam parameters_t The parameters type
  ///
  /// @param parameters TrackParameters to be checked
  /// @param bcheck BoundaryCheck directive for this onSurface check
  /// @return boolean indication if operation was successful
  template <typename parameters_t>
  bool
  onSurface(const parameters_t&  parameters,
            const BoundaryCheck& bcheck = BoundaryCheck(true)) const;

  /// The geometric onSurface method
  /// Geometrical check whether position is on Surface
  ///
  /// @param gpos global position to be evaludated
  /// @param bcheck BoundaryCheck directive for this onSurface check
  /// @return boolean indication if operation was successful
  virtual bool
  isOnSurface(const Vector3D& gpos, const BoundaryCheck& bcheck = true) const;

  /// The insideBounds method for local positions
  ///
  /// @param lpos local position to check
  /// @param bcheck  BoundaryCheck directive for this onSurface check
  /// @return boolean indication if operation was successful
  virtual bool
  insideBounds(const Vector2D& lpos, const BoundaryCheck& bcheck = true) const;

  /// Local to global transformation
  /// Generalized local to global transformation for the surface types. Since
  /// some surface types need the global momentum/direction to resolve sign
  /// ambiguity this is also provided
  ///
  /// @param lpos local 2D posittion in specialized surface frame
  /// @param gmom global 3D momentum representation (optionally ignored)
  /// @param gpos global 3D position to be filled (given by reference for method
  /// symmetry)
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& gmom,
                Vector3D&       gpos) const = 0;

  /// Global to local transformation
  /// Generalized global to local transformation for the surface types. Since
  /// some surface types need the global momentum/direction to resolve sign
  /// ambiguity this is also provided
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param gmom global 3D momentum representation (optionally ignored)
  /// @param lpos local 2D position to be filled (given by reference for method
  /// symmetry)
  /// @return boolean indication if operation was successful (fail means global
  /// position was not on surface)
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& gmom,
                Vector2D&       lpos) const = 0;

  /// Return mehtod for the reference frame
  /// This is the frame in which the covariance matrix is defined (specialized
  /// by all surfaces)
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param gmom global 3D momentum representation (optionally ignored)
  /// @return RotationMatrix3D which defines the three axes of the measurement
  /// frame
  virtual const Acts::RotationMatrix3D
  referenceFrame(const Vector3D& gpos, const Vector3D& gmom) const;

  /// Initialize the jacobian from local to global
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param jac is the jacobian to be initialized
  /// @param pos is the global position of the parameters
  /// @param dir is the direction at of the parameters
  /// @param pars is the parameter vector
  virtual void initJacobianToGlobal(ActsMatrixD<7, 5>& jac,
                                    const Vector3D&       gpos,
                                    const Vector3D&       dir,
                                    const ActsVectorD<5>& pars) const;

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
  virtual const RotationMatrix3D initJacobianToLocal(ActsMatrixD<5, 7>& jac,
                                                     const Vector3D& gpos,
                                                     const Vector3D& dir) const;

  /// Calculate the form factors for the derivatives
  /// the calculation is identical for all surfaces where the
  /// reference frame does not depend on the direction
  ///
  /// @param gpos is the position of the paramters in global
  /// @param dir is the direction of the track
  /// @param rft is the transposed reference frame (avoids recalculation)
  /// @param jac is the transport jacobian
  ///
  /// @return a five-dim vector
  virtual const ActsRowVectorD<5>
  derivativeFactors(const Vector3D&         gpos,
                    const Vector3D&         dir,
                    const RotationMatrix3D& rft,
                    const ActsMatrixD<7, 5>& jac) const;

  /// Calucation of the path correction for incident
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param gmom global 3D momentum representation
  /// @return Path correction with respect to the nominal incident.
  virtual double
  pathCorrection(const Vector3D& gpos, const Vector3D& gmom) const = 0;

  /// Straight line intersection schema from parameters
  ///
  /// Templated for :
  /// @tparam parameters_t Type of track parameters
  /// @tparam options_t Type of the navigation options
  /// @tparam corrector_t is the type of the corrector struct foer the direction
  ///
  /// @param parameters The parameters to start from
  /// @param options Options object that holds additional navigation info
  /// @param correct Corrector struct that can be used to refine the solution
  ///
  /// @return SurfaceIntersection object (contains intersection & surface)
  template <typename parameters_t,
            typename options_t,
            typename corrector_t = VoidCorrector>
  SurfaceIntersection
  intersectionEstimate(const parameters_t& parameters,
                       const options_t&    options,
                       const corrector_t&  correct = corrector_t()) const
  {
    // get the intersection with the surface
    auto sIntersection = intersectionEstimate(parameters.position(),
                                              parameters.direction(),
                                              options.navDir,
                                              options.boundaryCheck,
                                              correct);
    // return a surface intersection with result direction
    return SurfaceIntersection(sIntersection, this);
  }

  /// Straight line intersection from position and momentum
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  ///        inside bounds (check is done)
  /// @param 3D direction representation - expected to be normalized (no check
  /// done)
  /// @param navDir The navigation direction : if you want to find the closest,
  ///        chose anyDirection and the closest will be chosen
  /// @param bcheck boundary check directive for this operation
  /// @param corr is a correction function on position and momentum to do
  ///        a more appropriate intersection
  ///
  /// @return Intersection object
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gidr,
                       NavigationDirection  navDir = forward,
                       const BoundaryCheck& bcheck = false,
                       CorrFnc              corr = nullptr) const = 0;
  /// clang-format on

  /// Returns 'true' if this surface is 'free'
  /// i.e. it does not belong to a detector element, a layer or a tracking
  /// volume
  virtual bool
  isFree() const;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the ostream to be dumped into
  virtual std::ostream&
  dump(std::ostream& sl) const;

  /// Return properly formatted class name
  virtual std::string
  name() const = 0;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  virtual variant_data
  toVariantData() const = 0;

protected:
  /// Transform3D definition that positions
  /// (translation, rotation) the surface in global space
  std::shared_ptr<const Transform3D> m_transform;

  /// Pointer to the a DetectorElementBase
  const DetectorElementBase* m_associatedDetElement;

  /// Associated Identifier to this
  Identifier m_associatedDetElementId;

  /// The associated layer Layer - layer in which the Surface is be embedded,
  /// nullptr if not associated
  const Layer* m_associatedLayer;

  /// The assoicated TrackingVolume - tracking volume in case the surface is a
  /// boundary surface, nullptr if not
  const TrackingVolume* m_associatedTrackingVolume;

  /// Possibility to attach a material descrption
  std::shared_ptr<const SurfaceMaterial> m_associatedMaterial;
};

#include "detail/Surface.ipp"

std::ostream&
operator<<(std::ostream& sl, const Surface& sf);

}  // namespace Acts