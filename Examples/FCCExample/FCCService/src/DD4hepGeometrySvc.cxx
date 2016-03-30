#include "FCCService/DD4hepGeometrySvc.h"

DECLARE_SERVICE_FACTORY(Acts::DD4hepGeometrySvc)

Acts::DD4hepGeometrySvc::DD4hepGeometrySvc(const std::string& name, ISvcLocator* svc) :
Acts::ServiceBase(name, svc),
m_worldDetElement(),
m_log(msgSvc(), name)
{
    declareProperty("Detector", m_xmlFileNames, "XML file with detector description");
}

Acts::DD4hepGeometrySvc::~DD4hepGeometrySvc()
{}

StatusCode Acts::DD4hepGeometrySvc::initialize()
{
    MSG_DEBUG("initialize()");
    //Service needs to be initialized
    if (!ServiceBase::initialize()) return StatusCode::FAILURE;
    // retrieve the the static instance of the DD4HEP::Geometry
    DD4hep::Geometry::LCDD* lcdd = &(DD4hep::Geometry::LCDD::getInstance());
    lcdd->addExtension<IDD4hepGeometrySvc>(this);
    //load geometry from file
    for (auto& filename : m_xmlFileNames) {
        MSG_DEBUG("Loading DD4hep detector geometry from file:  '" << filename << "'");
        lcdd->fromCompact(filename);
    }
    lcdd->volumeManager();
    lcdd->apply("DD4hepVolumeManager",0,0);
    //set the world detector element
    m_worldDetElement = lcdd->world();
    
    return StatusCode::SUCCESS;
}

StatusCode Acts::DD4hepGeometrySvc::finalize() {
    return StatusCode::SUCCESS;
}

/** Query the interfaces. */
StatusCode Acts::DD4hepGeometrySvc::queryInterface(const InterfaceID& riid, void** ppvInterface)
{
    if ( Acts::IID_IDD4hepGeometrySvc == riid )
        *ppvInterface = (Acts::IDD4hepGeometrySvc*)this;
    else  {
        // Interface is not directly available: try out a base class
        return Service::queryInterface(riid, ppvInterface);
    }
    addRef();
    return StatusCode::SUCCESS;
}