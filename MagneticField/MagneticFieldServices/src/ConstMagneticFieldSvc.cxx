#include "MagneticFieldServices/ConstMagneticFieldSvc.h"

DECLARE_SERVICE_FACTORY(Acts::ConstMagneticFieldSvc)

Acts::ConstMagneticFieldSvc::ConstMagneticFieldSvc(const std::string& name, ISvcLocator* svc) :
Acts::ServiceBase(name,svc),
m_magField(0.)
{
    MSG_INFO("ConstMagneticFieldSvc Constructor");
    declareProperty("MagneticFieldValue",m_magField);
}

Acts::ConstMagneticFieldSvc::~ConstMagneticFieldSvc()
{}

StatusCode Acts::ConstMagneticFieldSvc::initialize()
{
    MSG_INFO("initilaize() ConstMagneticFieldSvc");
    //Service needs to be initialized
    if (!ServiceBase::initialize()) return StatusCode::FAILURE;
    return StatusCode::SUCCESS;
}

StatusCode Acts::ConstMagneticFieldSvc::finalize()
{
    MSG_DEBUG("finalize()");
    return StatusCode::SUCCESS;
}

StatusCode Acts::ConstMagneticFieldSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
    if ( IID_IMagneticFieldSvc == riid )
        *ppvInterface = (IMagneticFieldSvc*)this;
    else  {
        // Interface is not directly available: try out a base class
        return Service::queryInterface(riid, ppvInterface);
    }
    addRef();
    return StatusCode::SUCCESS;
}

void Acts::ConstMagneticFieldSvc::getField(const double *, double* bxyz,double *deriv) {
    //return on every point the same field value
    if (bxyz) for (int i=0; i<3; i++) bxyz[i]  = m_magField;
    else MSG_ERROR("Can not fill magnetic field, magnetic field pointer is nullptr");
    //if derivatives are asked return them
    if (deriv) for (int i=0; i<9; i++) deriv[i] = 0.;
}

void Acts::ConstMagneticFieldSvc::getFieldZR(const double *, double* bxyz,double *deriv) {
    //return on every point the same field value
    if (bxyz) for (int i=0; i<3; i++) bxyz[i]  = m_magField;
    else MSG_ERROR("Can not fill magnetic field, magnetic field pointer is nullptr");
    //if derivatives are asked return them
    if (deriv) for (int i=0; i<9; i++) deriv[i] = 0.;
}