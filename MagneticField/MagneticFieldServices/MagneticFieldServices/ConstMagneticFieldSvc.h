///////////////////////////////////////////////////////////////////
// ConstMagneticFieldSvc.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MAGNETICFIELDSERVICES_CONSTMAGFIELDSVC_H
#define ACTS_MAGNETICFIELDSERVICES_CONSTMAGFIELDSVC_H 1

// Core module
#include "CoreInterfaces/ServiceBase.h"
// Interface
#include "MagneticFieldInterfaces/IMagneticFieldSvc.h"

namespace Acts {
    
    /** @ class ConstMagneticFieldSvc
     
     Returns a const magnetic Field at every position. The cont value can be set over the job option file.
     
     @ author julia.hrdinka@cern.ch
     */
    class ConstMagneticFieldSvc: public Acts::ServiceBase, virtual public IMagneticFieldSvc {
        
    public:
        /** constructor */
        ConstMagneticFieldSvc(const std::string& name, ISvcLocator* svc);
        /** destructor */
        ~ConstMagneticFieldSvc();
        /** Framework methods */
        virtual StatusCode initialize() final;
        virtual StatusCode finalize() final;
        
        /** Query the interfaces.
         /   Input: riid, Requested interface ID
         /          ppvInterface, Pointer to requested interface **/
        StatusCode queryInterface( const InterfaceID& riid, void** ppvInterface );
        
        /** Retrieve interface ID */
        static const InterfaceID& interfaceID() { return IID_IMagneticFieldSvc;}
        
        /** get B field value at given position */
        /** xyz[3] is in mm, bxyz[3] is in kT */
        /** if deriv[9] is given, field derivatives are returned in kT/mm */
        virtual void getField(const double *xyz, double *bxyz, double *deriv = nullptr) final;
        
        /** get B field value on the z-r plane at given position */
        /** works only inside the solenoid; otherwise calls getField() above */
        /** xyz[3] is in mm, bxyz[3] is in kT */
        /** if deriv[9] is given, field derivatives are returned in kT/mm */
        virtual void getFieldZR(const double *xyz, double *bxyz, double *deriv = nullptr) final;
        
    private:
        
        double                          m_magField; //!<value of the magnetic field
        
    };
}


#endif //>  ACTS_MAGNETICFIELDSERVICES_CONSTMAGFIELDSVC_H

