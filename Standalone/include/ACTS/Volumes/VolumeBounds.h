///////////////////////////////////////////////////////////////////
// VolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_VOLUMEBOUNDS_H
#define ACTS_VOLUMES_VOLUMEBOUNDS_H 1

// STD/STL
#include <iostream>
#include <iomanip>
#include <memory>

// Geometry module
#include "ACTS/GeometryUtils/BinningType.h"
// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"

namespace Acts {


  class Surface;
  class Volume;

  // master typedef
  class VolumeBounds;
  typedef std::shared_ptr<const VolumeBounds> VolumeBoundsPtr;

  /**
   @class VolumeBounds

   Pure Absract Base Class for Volume bounds.

   Acts::VolumeBounds are a set of up to six confining Surfaces that are stored in a std::vector.
   Each type of Acts::VolumeBounds has to implement a decomposeToSurfaces() and a inside() method.

   The orientation of the Surfaces are in a way that the normal vector points to the oustide world.

   The Volume, retrieving a set of Surfaces from the VolumeBounds, can turn the Surfaces into BoundarySurfaces.

   @author Andreas.Salzburger@cern.ch
   */

  class VolumeBounds {

    public:
      /**Default Constructor*/
      VolumeBounds(){}

      /**Destructor*/
      virtual ~VolumeBounds(){}

      /** clone() method to make deep copy in Volume copy constructor and for assigment operator
         of the Surface class.*/
      virtual VolumeBounds* clone() const = 0;

      /** Checking if position given in volume frame is inside*/
      virtual bool inside(const Vector3D& pos, double tol=0.) const = 0;

      /** Method to decompose the Bounds into Surfaces, the Volume can turn them into BoundarySurfaces */
      virtual const std::vector<const Acts::Surface*>* decomposeToSurfaces(std::shared_ptr<Transform3D> transform) const = 0;

      /** Binning offset - overloaded for some R-binning types */
      virtual Vector3D binningOffset(BinningValue bValue) const;

      /** Binning borders in double */
      virtual double binningBorder(BinningValue bValue) const;

      /** Output Method for std::ostream, to be overloaded by child classes */
      virtual std::ostream& dump(std::ostream& sl) const = 0;

  };

  /** Binning offset - overloaded for some R-binning types */
  inline Vector3D VolumeBounds::binningOffset(BinningValue) const
  { // standard offset is 0.,0.,0.
    return Vector3D(0.,0.,0.);
  }

  inline double VolumeBounds::binningBorder(BinningValue) const
  {
     return 0.;
  }

  /**Overload of << operator for std::ostream for debug output*/
  std::ostream& operator << ( std::ostream& sl, const VolumeBounds& vb);

} // end of namespace

#endif // ACTS_VOLUMES_VOLUMEBOUNDS_H
