///////////////////////////////////////////////////////////////////
// MagneticFieldProperties.cxx, ATS project
///////////////////////////////////////////////////////////////////

// MagneticField module
#include "MagneticFieldUtils/MagneticFieldProperties.h"

Ats::MagneticFieldProperties::MagneticFieldProperties(Ats::MagneticFieldMode mode) :
  m_magneticFieldMode(mode),
  m_magneticField(0., 0., 0.)
{}

Ats::MagneticFieldProperties::MagneticFieldProperties(const Vector3D& field) :
  m_magneticFieldMode(Ats::ConstantField),
  m_magneticField(field)
{}


Ats::MagneticFieldProperties::MagneticFieldProperties(const Ats::MagneticFieldProperties& magprop) :
  m_magneticFieldMode(magprop.m_magneticFieldMode),
  m_magneticField(magprop.m_magneticField)
{}  

Ats::MagneticFieldProperties& Ats::MagneticFieldProperties::operator=(const Ats::MagneticFieldProperties& magprop)
{
  if (this != &magprop){
     m_magneticFieldMode = magprop.m_magneticFieldMode;
     m_magneticField     = magprop.m_magneticField;
  }
  return(*this);
}
 
/**Overload of << operator for both, MsgStream and std::ostream for debug output*/ 
MsgStream& Ats::operator<<( MsgStream& sl, const Ats::MagneticFieldProperties& mprop)
{
   sl << "Ats::MagneticFieldProperties, configuration: " << mprop.magneticFieldMode() << endreq; 
   return sl;
}

std::ostream& Ats::operator << ( std::ostream& sl, const Ats::MagneticFieldProperties& mprop)
{
   sl << "Ats::MagneticFieldProperties, configuration: " << mprop.magneticFieldMode() << std::endl;
   return sl;
}
