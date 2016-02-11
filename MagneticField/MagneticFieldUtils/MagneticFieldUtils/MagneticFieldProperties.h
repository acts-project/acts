///////////////////////////////////////////////////////////////////
// MagneticFieldProperties.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_MAGNETICFIELDUTILS_MAGNETICFIELDPROPERTIES_H
#define ATS_MAGNETICFIELDUTILS_MAGNETICFIELDPROPERTIES_H 1

// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/GaudiException.h"
// MagneticField module
#include "MagneticFieldUtils/MagneticFieldMode.h"
// Core module
#include "Algebra/AlgebraDefinitions.h"
// STD
#include <iostream>

namespace Ats {

  /** 
   @class MagneticFieldProperties

     new magnetic field properties to steer the behavior of the extrapolation
     
   @author Andreas.Salzburger@cern.ch 
   */
  class MagneticFieldProperties {
      
    public:
      /**Constructor for magnetic field mode - full field is default */
      MagneticFieldProperties(MagneticFieldMode mode=Ats::FullField);

      /**Constructor for magnetic field mode */
      MagneticFieldProperties(const Vector3D& field);
      
      /**Copy Constructor */  
      MagneticFieldProperties(const MagneticFieldProperties& matprop);
    
      /**Destructor*/
      virtual ~MagneticFieldProperties(){}
      
      /**Assignment operator */
      MagneticFieldProperties& operator=(const MagneticFieldProperties& matprop);
      
      /**Cast operator*/
      operator MagneticFieldMode () const;
            
      /**Returns the MagneticFieldMode as specified */
      MagneticFieldMode magneticFieldMode() const;
            
      /** Get the magnetic field - in case of constant field only - throws exception if mode is not constant */
      const Vector3D& magneticField() const throw (GaudiException);

    protected:
      mutable MagneticFieldMode         m_magneticFieldMode;
      Vector3D                          m_magneticField;
  };


  inline MagneticFieldProperties::operator MagneticFieldMode() const { return m_magneticFieldMode; }  

  inline MagneticFieldMode MagneticFieldProperties::magneticFieldMode() const { return m_magneticFieldMode; }  
  
  inline const Vector3D& MagneticFieldProperties::magneticField() const throw (GaudiException) { 
      if ( m_magneticFieldMode != Ats::ConstantField ) 
          throw GaudiException("Ats::MagneticFieldProperties", "You can only ask for a field value if you have a constant field!", StatusCode::FAILURE);
      return m_magneticField;
  }


/**Overload of << operator for both, MsgStream and std::ostream for debug output*/ 
MsgStream& operator << ( MsgStream& sl, const MagneticFieldProperties& mprop);

std::ostream& operator << ( std::ostream& sl, const MagneticFieldProperties& mprop);
    
} // end of namespace

#endif // TRKGEOMETRY_MAGNETICFIELDPROPERTIES_H


