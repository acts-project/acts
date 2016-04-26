///////////////////////////////////////////////////////////////////
// MagneticFieldProperties.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MAGNETICFIELDUTILS_MAGNETICFIELDPROPERTIES_H
#define ACTS_MAGNETICFIELDUTILS_MAGNETICFIELDPROPERTIES_H 1

#include "ACTS/MagneticField/MagneticFieldMode.h"
#include "ACTS/Utilities/Definitions.h"
#include <iostream>

namespace Acts {

  /** 
   @class MagneticFieldProperties

     new magnetic field properties to steer the behavior of the extrapolation
     
   @author Andreas.Salzburger@cern.ch 
   */
  class MagneticFieldProperties {
      
    public:
      /**Constructor for magnetic field mode - full field is default */
      MagneticFieldProperties(MagneticFieldMode mode=Acts::FullField);

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
      const Vector3D& magneticField() const;

    protected:
      mutable MagneticFieldMode         m_magneticFieldMode;
      Vector3D                          m_magneticField;
  };


  inline MagneticFieldProperties::operator MagneticFieldMode() const { return m_magneticFieldMode; }  

  inline MagneticFieldMode MagneticFieldProperties::magneticFieldMode() const { return m_magneticFieldMode; }  
  
  inline const Vector3D& MagneticFieldProperties::magneticField() const { 
      return m_magneticField;
  }


/**Overload of << operator for std::ostream */ 
  std::ostream& operator << ( std::ostream& sl, const MagneticFieldProperties& mprop);
    
} // end of namespace

#endif // ACTS_MAGNETICFIELDUTILS_MAGNETICFIELDPROPERTIES_H


