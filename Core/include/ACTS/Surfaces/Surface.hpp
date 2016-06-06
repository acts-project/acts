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

/**
 @class Surface

 Abstract Base Class for tracking surfaces

 The creation of a Surface by passing a std::shared_ptr<Transform3D> through the
 constructor
 implies that the ownership of the Transform3D object is also passed to the
 Surface,
 therfor the memory is freed in the Surface destructor.

 For all isOnSurface, or positionOnSurface and insideBounds methods two
 tolerance parameters
 can be given which correspond to the two local natural coordinates of the
 surface loc1, loc2.

 */

class Surface : public virtual GeometryObject
{
public:
  /** @enum SurfaceType

      This enumerator simplifies the persistency & calculations,
      by saving a dynamic_cast to happen.

      Other is reserved for the GeometrySurfaces implementation.

    */
  enum SurfaceType {
    Cone        = 0,
    Cylinder    = 1,
    Disc        = 2,
    Perigee     = 3,
    Plane       = 4,
    Line        = 5,
    Curvilinear = 6,
    Other       = 7
  };

  /** Default Constructor - needed for inherited classes */
  Surface();

  /** Copy constructor - it invalidates the association to detector element and
   * identifier */
  Surface(const Surface& sf);

  /** Copy constructor with shift */
  Surface(const Surface& sf, const Transform3D& transf);

  /** Constructor with Transform3D, surface often share transforms */
  Surface(std::shared_ptr<Transform3D> htrans);

  /** Constructor with Transform3D, by unique_ptr */
  Surface(std::unique_ptr<Transform3D> htrans);

  /** Constructor from TrkDetElement*/
  Surface(const DetectorElementBase& detelement);

  /** Constructor form TrkDetElement and Identifier*/
  Surface(const DetectorElementBase& detelement, const Identifier& id);

  /** Virtual Destructor*/
  virtual ~Surface();

  /** Assignment operator - the alssignment invalidates the link
     to detector element and identifier */
  Surface&
  operator=(const Surface& sf);

  /** Equality operator*/
  virtual bool
  operator==(const Surface& sf) const = 0;

  /** Non-equality operator*/
  virtual bool
  operator!=(const Surface& sf) const;

  /** Implicit constructor - uses the copy constructor
    - a new position can optionally be given  */
  virtual Surface*
  clone(const Transform3D* shift = nullptr) const = 0;

  /** Returns the Surface type to avoid dynamic casts */
  virtual SurfaceType
  type() const = 0;

  /** Return the cached transformation directly, will create
      a new transform if it's not here. */
  std::shared_ptr<Transform3D>
  cachedTransform() const;

  /** Returns Transform3D by reference */
  virtual const Transform3D&
  transform() const;

  /** Returns the center position of the Surface */
  virtual const Vector3D&
  center() const;

  /** Returns the normal vector of the Surface */
  virtual const Vector3D&
  normal() const;

  /** Returns a normal vector at a specific local position */
  virtual const Vector3D
  normal(const Vector2D& lp) const;

  /** Returns a normal vector at a specific local position - no checking done on
   * global position */
  virtual const Vector3D
  normal(const Vector3D& global) const;

  /** The binning position method - as default the center is given, but may be
   * overloaded */
  virtual Vector3D
  binningPosition(BinningValue bValue) const override;

  /** return associated Detector Element */
  const DetectorElementBase*
  associatedDetectorElement() const;

  /** return Identifier of the associated Detector Element */
  const Identifier
  associatedDetectorElementIdentifier() const;

  /** return the associated Layer */
  const Layer*
  associatedLayer() const;

  /** The templated Parameters OnSurface method - checks on surface pointer
   * first */
  template <class T>
  bool
  onSurface(const T&             parameters,
            const BoundaryCheck& bchk = BoundaryCheck(true)) const;

  /** This method returns true if the GlobalPosition is on the Surface for both,
    within
    or without check of whether the local position is inside boundaries or not
    */
  virtual bool
  isOnSurface(const Vector3D& glopo, const BoundaryCheck& bchk = true) const;

  /**  virtual methods to be overwritten by the inherited surfaces */
  virtual bool
  insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk) const;

  /** Specified by each surface type: LocalToGlobal method without dynamic
   * memory allocation */
  virtual void
  localToGlobal(const Vector2D& locp,
                const Vector3D& mom,
                Vector3D&       glob) const = 0;

  /** Specified by each surface type: GlobalToLocal method without dynamic
   * memory allocation - boolean checks if on surface */
  virtual bool
  globalToLocal(const Vector3D& glob,
                const Vector3D& mom,
                Vector2D&       loc) const = 0;

  /** the pathCorrection for derived classes with thickness - it reflects if the
   * direction projection is positive or negative */
  virtual double
  pathCorrection(const Vector3D& pos, const Vector3D& mom) const;

  /** Return the measurement frame - this is needed for alignment, in particular
     for StraightLine and Perigee Surface
       - the default implementation is the the RotationMatrix3D of the transform
     */
  virtual const Acts::RotationMatrix3D
  measurementFrame(const Vector3D& glopos, const Vector3D& glomom) const;

  /** fst straight line intersection schema - templated for cvharged and neutral
   * parameters */
  template <class T>
  Intersection
  intersectionEstimate(const T&            pars,
                       bool                forceDir = false,
                       Acts::BoundaryCheck bchk     = false) const
  {
    return intersectionEstimate(
        pars.position(), pars.momentum().unit(), forceDir, bchk);
  }

  /** fast straight line intersection schema - standard: provides closest
     intersection and (signed) path length
      forceFwd is to provide the closest forward solution */
  virtual Intersection
  intersectionEstimate(const Vector3D&      pos,
                       const Vector3D&      dir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bchk = false) const = 0;

  /** Surface Bounds method */
  virtual const SurfaceBounds&
  bounds() const = 0;

  /** Returns 'true' if this surface is 'free', i.e. it does not belong to a
   * detector element (and returns false otherwise*/
  bool
  isFree() const;

  void
  setSurfaceMaterial(std::shared_ptr<const SurfaceMaterial> material) const;

  /** Return 'true' if this surface is own by the detector element */
  bool
  isActive() const;

  /** set material layer */
  const SurfaceMaterial*
  surfaceMaterial() const;

  /** Output Method for std::ostream, to be overloaded by child classes */
  virtual std::ostream&
  dump(std::ostream& sl) const;

  /** Return properly formatted class name */
  virtual std::string
  name() const = 0;

  /**return number of surfaces currently created - needed for EDM monitor */
  static unsigned int
  numberOfInstantiations();

  /**return number of free surfaces currently created (i.e. those not belonging
   * to a DE) - needed for EDM monitor */
  static unsigned int
  numberOfFreeInstantiations();

  /** method to associate the associated Acts::Layer which is alreay owned
     - only allowed by LayerBuilder
     - only done if no Layer is set already  */
  void
  associateLayer(const Layer&) const;

protected:
  /** Private members are in principle implemented as mutable pointers to
    objects for easy checks
    if they are already declared or not */
  std::shared_ptr<Transform3D>
                    m_transform;  //!< Transform3D w.r.t to global frame
  mutable Vector3D* m_center;     //!< center position of the surface
  mutable Vector3D* m_normal;     //!< normal vector of the surface

  /** Pointer to the a DetectorElementBase */
  const DetectorElementBase* m_associatedDetElement;
  Identifier                 m_associatedDetElementId;

  /** The associated layer Acts::Layer
   - layer in which the Surface is be embedded */
  mutable const Layer* m_associatedLayer;

  /** Possibility to attach a material descrption */
  mutable std::shared_ptr<const SurfaceMaterial> m_surfaceMaterial;

  /** number of objects of this type in memory - needed for EDM monitor*/
  static unsigned int s_numberOfInstantiations;

  /** number of objects of this type in memory which do not belong to a detector
   * element - needed for EDM monitor*/
  static unsigned int s_numberOfFreeInstantiations;
};

inline bool
Surface::operator!=(const Surface& sf) const
{
  return !((*this) == sf);
}

inline std::shared_ptr<Transform3D>
Surface::cachedTransform() const
{
  return m_transform;
}

inline const Transform3D&
Surface::transform() const
{
  if (m_transform.get()) return (*(m_transform.get()));
  if (m_associatedDetElement && m_associatedDetElementId.is_valid())
    return m_associatedDetElement->transform(m_associatedDetElementId);
  if (m_associatedDetElement) return m_associatedDetElement->transform();
  return s_idTransform;
}

inline const Vector3D&
Surface::center() const
{
  if (m_transform.get() && !m_center)
    m_center = new Vector3D(m_transform->translation());
  if (m_center) return (*m_center);
  if (m_associatedDetElement && m_associatedDetElementId.is_valid())
    return m_associatedDetElement->center(m_associatedDetElementId);
  if (m_associatedDetElement) return m_associatedDetElement->center();
  return s_origin;
}

inline const Vector3D&
Surface::normal() const
{
  if (m_transform.get() && m_normal == 0)
    m_normal = new Vector3D(m_transform->rotation().col(2));
  if (m_normal) return (*m_normal);
  if (m_associatedDetElement && m_associatedDetElementId.is_valid())
    return m_associatedDetElement->normal(m_associatedDetElementId);
  if (m_associatedDetElement) return m_associatedDetElement->normal();
  return s_zAxis;
}

// return the binning position for ordering in the BinnedArray
inline Vector3D Surface::binningPosition(BinningValue) const
{
  // very simple binning directives following hte binning type
  // give the center as default for all of these binning types
  // binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
  return center();
}

// common to planar surfaces
inline double
Surface::pathCorrection(const Vector3D&, const Vector3D& mom) const
{
  Vector3D dir(mom.unit());
  double   cosAlpha = dir.dot(normal());
  return fabs(1. / cosAlpha);
}

//* the templated parameters on Surface method */
template <class T>
bool
Surface::onSurface(const T& pars, const BoundaryCheck& bcheck) const
{
  // surface pointer comparison as a first fast check (w/o transform)
  if ((&pars.associatedSurface() == this)) {
    return (bcheck ? insideBounds(pars.localPosition(), bcheck) : true);
  }
  return isOnSurface(pars.position(), bcheck);
}

inline bool
Surface::insideBounds(const Vector2D& locpos, const BoundaryCheck& bchk) const
{
  return bounds().inside(locpos, bchk);
}

// take local position and return a normal direction, local position is ignored
// for planar surfaces
inline const Vector3D
Surface::normal(const Vector2D&) const
{
  return normal();
}

// take local position and return a normal direction, local position is ignored
// for planar surfaces
inline const Vector3D
Surface::normal(const Vector3D&) const
{
  return normal();
}

inline const DetectorElementBase*
Surface::associatedDetectorElement() const
{
  return m_associatedDetElement;
}

inline const Identifier
Surface::associatedDetectorElementIdentifier() const
{
  if (!m_associatedDetElement) return Identifier();  // in invalid state
  if (m_associatedDetElementId.is_valid()) return m_associatedDetElementId;
  return m_associatedDetElement->identify();
}

inline const Layer*
Surface::associatedLayer() const
{
  return (m_associatedLayer);
}

inline const SurfaceMaterial*
Surface::surfaceMaterial() const
{
  return m_surfaceMaterial.get();
}

inline bool
Surface::isActive() const
{
  return (m_associatedDetElement != nullptr);
}

inline bool
Surface::isFree() const
{
  return (m_associatedDetElement == nullptr && m_associatedLayer == nullptr);
}

inline void
Surface::setSurfaceMaterial(
    std::shared_ptr<const SurfaceMaterial> material) const
{
  m_surfaceMaterial = material;
}

inline void
Surface::associateLayer(const Layer& lay) const
{
  m_associatedLayer = (&lay);
}

std::ostream&
operator<<(std::ostream& sl, const Surface& sf);

/** Surfaces are not constructed by shared_ptr factories as they are EDM
   objects,
    hence the SurfaceArray is defined as such: */
typedef BinnedArray<const Surface*> SurfaceArray;

}  // end of namespace Acts

#endif  // ACTS_SURFACES_SURFACE_H
