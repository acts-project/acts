// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Surface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_SURFACE_H
#define ACTS_SURFACES_SURFACE_H 1

#include "ACTS/Detector/DetectorElementBase.hpp"
#include "ACTS/Surfaces/BoundaryCheck.hpp"
#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometryObject.hpp"
#include "ACTS/Utilities/GeometryStatics.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include "ACTS/Utilities/Intersection.hpp"

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
/// The Surface class builds the core of the ACTS Tracking Geometry.
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
  /// surface pointer first.
  /// If that check fails, it calls the geometrical check isOnSurface
  ///
  /// @param parameters TrackParameters to be checked
  /// @param bcheck BoundaryCheck directive for this onSurface check
  /// @return boolean indication if operation was successful
  template <class T>
  bool
  onSurface(const T&             parameters,
            const BoundaryCheck& bcheck = BoundaryCheck(true)) const;

  /// The insideBounds method for local positions
  ///
  /// @param lpos local position to check
  /// @param bcheck  BoundaryCheck directive for this onSurface check
  /// @return boolean indication if operation was successful
  virtual bool
  insideBounds(const Vector2D& lpos, const BoundaryCheck& bcheck = true) const;

  /// The geometric onSurface method
  /// Geometrical check whether position is on Surface
  ///
  /// @param gpos global position to be evaludated
  /// @param bcheck BoundaryCheck directive for this onSurface check
  /// @return boolean indication if operation was successful
  virtual bool
  isOnSurface(const Vector3D& gpos, const BoundaryCheck& bcheck = true) const;

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

  /// Calucation of the path correction for incident
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param gmom global 3D momentum representation
  /// @return Path correction with respect to the nominal incident.
  virtual double
  pathCorrection(const Vector3D& gpos, const Vector3D& gmom) const = 0;

  /// Straight line intersection schema from parameters
  /// Templated for charged and neutral
  ///
  /// @todo include intersector
  /// @param pars TrackParameters to start from
  /// @param forceDir boolean indication whether to force the direction given by
  /// the TrackParameters to hold
  /// @param bcheck boundary check directive for this operation
  /// @return Intersection class
  template <class T>
  Intersection
  intersectionEstimate(const T&             pars,
                       bool                 forceDir = false,
                       const BoundaryCheck& bcheck   = false) const
  {
    return intersectionEstimate(
        pars.position(), pars.momentum().unit(), forceDir, bcheck);
  }

  /// Straight line intersection schema from parameters
  /// Templated for charged and neutral
  /// @todo include intersector
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  ///         inside bounds (check is done)
  /// @param gdir global 3D direction representation
  ///        @note has to be normalized
  /// @param forceDir boolean indication whether to force the direction given by
  /// the TrackParameters to hold
  /// @param bcheck boundary check directive for this operation
  ///
  /// @return Intersection class
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gdir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bcheck = false) const = 0;

  /// Distance to intersection
  template <class T>
  double
  distance(const T&             pars,
           bool                 forceDir = false,
           const BoundaryCheck& bcheck   = false) const;

  /// Returns 'true' if this surface is 'free'
  /// i.e. it does not belong to a detector element, a layer or a tracking
  /// volume
  bool
  isFree() const;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the ostream to be dumped into
  virtual std::ostream&
  dump(std::ostream& sl) const;

  /// Return properly formatted class name
  virtual std::string
  name() const = 0;

protected:
  /// Transform3D definition that positions (translation, rotation) the surface
  /// in global space
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

inline const RotationMatrix3D
Surface::referenceFrame(const Acts::Vector3D&, const Acts::Vector3D&) const
{
  return transform().matrix().block<3, 3>(0, 0);
}

template <class T>
bool
Surface::onSurface(const T& pars, const BoundaryCheck& bcheck) const
{
  // surface pointer comparison as a first fast check (w/o transform)
  if ((&pars.referenceSurface() == this)) {
    return (bcheck ? insideBounds(pars.localPosition(), bcheck) : true);
  }
  return isOnSurface(pars.position(), bcheck);
}

inline const DetectorElementBase*
Surface::associatedDetectorElement() const
{
  return m_associatedDetElement;
}

inline const Identifier
Surface::associatedIdentifier() const
{
  if (!m_associatedDetElement) return Identifier();  // in invalid state
  if (m_associatedDetElementId.is_valid()) return m_associatedDetElementId;
  return m_associatedDetElement->identify();
}

inline bool
Surface::isFree() const
{
  return (!m_associatedDetElement && !m_associatedTrackingVolume);
}

inline const Layer*
Surface::associatedLayer() const
{
  return (m_associatedLayer);
}

inline const SurfaceMaterial*
Surface::associatedMaterial() const
{
  return m_associatedMaterial.get();
}

inline void
Surface::setAssociatedMaterial(std::shared_ptr<const SurfaceMaterial> material)
{
  m_associatedMaterial = material;
}

inline void
Surface::associateLayer(const Layer& lay)
{
  m_associatedLayer = (&lay);
}

std::ostream&
operator<<(std::ostream& sl, const Surface& sf);

}  // end of namespace Acts

#endif  // ACTS_SURFACES_SURFACE_H
