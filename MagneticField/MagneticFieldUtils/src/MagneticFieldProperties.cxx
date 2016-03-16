///////////////////////////////////////////////////////////////////
// MagneticFieldProperties.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// MagneticField module
#include "MagneticFieldUtils/MagneticFieldProperties.h"

Acts::MagneticFieldProperties::MagneticFieldProperties(Acts::MagneticFieldMode mode) :
  m_magneticFieldMode(mode),
  m_magneticField(0., 0., 0.)
{}

Acts::MagneticFieldProperties::MagneticFieldProperties(const Vector3D& field) :
  m_magneticFieldMode(Acts::ConstantField),
  m_magneticField(field)
{}


Acts::MagneticFieldProperties::MagneticFieldProperties(const Acts::MagneticFieldProperties& magprop) :
  m_magneticFieldMode(magprop.m_magneticFieldMode),
  m_magneticField(magprop.m_magneticField)
{}  

Acts::MagneticFieldProperties& Acts::MagneticFieldProperties::operator=(const Acts::MagneticFieldProperties& magprop)
{
  if (this != &magprop){
     m_magneticFieldMode = magprop.m_magneticFieldMode;
     m_magneticField     = magprop.m_magneticField;
  }
  return(*this);
}
 
/**Overload of << operator for both, MsgStream and std::ostream for debug output*/ 
MsgStream& Acts::operator<<( MsgStream& sl, const Acts::MagneticFieldProperties& mprop)
{
   sl << "Acts::MagneticFieldProperties, configuration: " << mprop.magneticFieldMode() << endreq; 
   return sl;
}

std::ostream& Acts::operator << ( std::ostream& sl, const Acts::MagneticFieldProperties& mprop)
{
   sl << "Acts::MagneticFieldProperties, configuration: " << mprop.magneticFieldMode() << std::endl;
   return sl;
}
