///////////////////////////////////////////////////////////////////
// Surface.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/Surface.h"

// STD/STL
#include <iostream>
#include <iomanip>

unsigned int Acts::Surface::s_numberOfInstantiations     = 0;
unsigned int Acts::Surface::s_numberOfFreeInstantiations = 0;

Acts::Surface::Surface() :
  m_transform(nullptr),
  m_center(nullptr),
  m_normal(nullptr),
  m_associatedDetElement(nullptr),
  m_associatedDetElementId(),
  m_associatedLayer(nullptr),
  m_surfaceMaterial(nullptr)
{
#ifndef NDEBUG
  s_numberOfInstantiations++;     // EDM Monitor
  s_numberOfFreeInstantiations++; // EDM Monitor
#endif
}

Acts::Surface::Surface(std::shared_ptr<Acts::Transform3D> tform) :
  m_transform(tform),
  m_center(nullptr),
  m_normal(nullptr),
  m_associatedDetElement(nullptr),
  m_associatedDetElementId(),
  m_associatedLayer(nullptr),
  m_surfaceMaterial(nullptr)
{
	m_transform = tform;
#ifndef NDEBUG
  s_numberOfInstantiations++; // EDM Monitor - increment one instance
  s_numberOfFreeInstantiations++;
#endif
}

Acts::Surface::Surface(std::unique_ptr<Acts::Transform3D> tform) :
  Surface (std::shared_ptr<Acts::Transform3D>(tform.release()))
{
  // No EDM monitor here since we delegate to the previous constructor.
}

Acts::Surface::Surface(const Acts::DetectorElementBase& detelement) :
  m_transform(nullptr),
  m_center(nullptr),
  m_normal(nullptr),
  m_associatedDetElement(&detelement),
  m_associatedDetElementId(),
  m_associatedLayer(nullptr),
  m_surfaceMaterial(nullptr)
{
#ifndef NDEBUG
  s_numberOfInstantiations++; // EDM Monitor - increment one instance
#endif
}

Acts::Surface::Surface(const Acts::DetectorElementBase& detelement, const Identifier& id) :
  m_transform(nullptr),
  m_center(nullptr),
  m_normal(nullptr),
  m_associatedDetElement(&detelement),
  m_associatedDetElementId(id),
  m_associatedLayer(nullptr),
  m_surfaceMaterial(nullptr)
{
#ifndef NDEBUG
  s_numberOfInstantiations++; // EDM Monitor - increment one instance
#endif
}

// copy constructor - Attention! sets the associatedDetElement to 0 and the identifier to invalid
Acts::Surface::Surface(const Surface& sf) :
  m_transform(sf.m_transform),
  m_center(nullptr),
  m_normal(nullptr),
  m_associatedDetElement(nullptr),
  m_associatedDetElementId(),
  m_associatedLayer(sf.m_associatedLayer),
  m_surfaceMaterial(sf.m_surfaceMaterial)
{
#ifndef NDEBUG
  s_numberOfInstantiations++; // EDM Monitor - increment one instance
  // this is by definition a free surface since a copy is not allowed to point to the det element
  s_numberOfFreeInstantiations++;
#endif
}

// copy constructor with shift - Attention! sets the associatedDetElement to 0 and the identifieer to invalid
// also invalidates the material layer
Acts::Surface::Surface(const Surface& sf, const Acts::Transform3D& shift) :
  m_transform(std::shared_ptr<Acts::Transform3D>(new Acts::Transform3D(shift*sf.transform()))),
  m_center(nullptr),
  m_normal(nullptr),
  m_associatedDetElement(nullptr),
  m_associatedDetElementId(),
  m_associatedLayer(nullptr),
  m_surfaceMaterial(nullptr)
{
#ifndef NDEBUG
  s_numberOfInstantiations++; // EDM Monitor - increment one instance
  // this is by definition a free surface since a copy is not allowed to point to the det element
  s_numberOfFreeInstantiations++;
#endif
}

// destructor
Acts::Surface::~Surface()
{
#ifndef NDEBUG
  s_numberOfInstantiations--; // EDM Monitor - decrement one instance
  if ( isFree() ) s_numberOfFreeInstantiations--;
#endif

  delete m_center;
  delete m_normal;
}

// assignment operator
// the assigned surfaces loses its link to the detector element
Acts::Surface& Acts::Surface::operator=(const Acts::Surface& sf)
{
  if (this!=&sf){
    delete m_center;  m_center = 0;
    delete m_normal;  m_normal = 0;
    m_transform              = sf.m_transform;
    m_associatedDetElement   = nullptr; // link to detector element is forced to be broken
    m_associatedDetElementId = Identifier();
    m_associatedLayer        = sf.m_associatedLayer;
    m_surfaceMaterial        = sf.m_surfaceMaterial;
  }
  return *this;
}

// checks if GlobalPosition is on Surface and inside bounds
bool Acts::Surface::isOnSurface(const Acts::Vector3D& glopo, const BoundaryCheck& bchk) const
{
    // create the local position
    Acts::Vector2D lpos;
    // global to local transformation
    bool g2L = globalToLocal(glopo,Acts::Vector3D::UnitX(),lpos);
    if (g2L) {
        // no boundary check, then return true
        if (!bchk) return true;
        // return what ever the bounds tell you
        return bounds().inside(lpos,bchk);
    }
    // did not succeed
    return false;
}

// return the measurement frame
const Acts::RotationMatrix3D Acts::Surface::measurementFrame(const Acts::Vector3D&, const Acts::Vector3D&) const
{
  return transform().rotation();
}

// for the EDM monitor
unsigned int Acts::Surface::numberOfInstantiations()
{
  return s_numberOfInstantiations;
}

// for the EDM monitor
unsigned int Acts::Surface::numberOfFreeInstantiations()
{
    return s_numberOfFreeInstantiations;
}

// overload dump for stream operator
std::ostream& Acts::Surface::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(4);
    sl << name() << std::endl;
    sl << "     Center position  (x, y, z) = (" << center().x() << ", " << center().y() << ", " << center().z() << ")" << std::endl;
    Acts::RotationMatrix3D rot(transform().rotation());
    Acts::Vector3D  rotX(rot.col(0));
    Acts::Vector3D  rotY(rot.col(1));
    Acts::Vector3D  rotZ(rot.col(2));
    sl << std::setprecision(6);
    sl << "     Rotation:             colX = (" << rotX(0) << ", " << rotX(1) << ", " << rotX(2) << ")" << std::endl;
    sl << "                           colY = (" << rotY(0) << ", " << rotY(1) << ", " << rotY(2) << ")" << std::endl;
    sl << "                           colZ = (" << rotZ(0) << ", " << rotZ(1) << ", " << rotZ(2) << ")" << std::endl;
    sl << "     Bounds  : " << bounds();
    sl << std::setprecision(-1);
    return sl;
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream& Acts::operator << ( std::ostream& sl, const Acts::Surface& sf)
{ return sf.dump(sl); }

