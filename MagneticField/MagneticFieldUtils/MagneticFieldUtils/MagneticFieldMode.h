///////////////////////////////////////////////////////////////////
// MagneticFieldMode.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_MAGNETICFIELDUTILS_MAGNETICFIELDMODE_H
#define ATS_MAGNETICFIELDUTILS_MAGNETICFIELDMODE_H 1

namespace Ats {

      /** @enum MagneticFieldMode

          MagneticFieldMode describing the field setup within a volume
        
         @author Andreas.Salzburger@cern.ch
        */
      enum  MagneticFieldMode { NoField                 = 0,  //!< Field is set to 0., 0., 0.,
                                ConstantField           = 1,  //!< Field is set to be constant
                                FastField               = 2,  //!< call the fast field access method of the FieldSvc
                                FullField               = 3   //!< Field is set to be realistic, but within a given Volume                                
                                };
      
} // end of namespace
#endif // ATS_MAGNETICFIELDUTILS_MAGNETICFIELDMODE_H
