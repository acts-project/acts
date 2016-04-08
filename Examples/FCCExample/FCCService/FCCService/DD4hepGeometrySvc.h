#ifndef FCCGEOMETRYSERVICES_DD4HEPGEOMETRYSVC_H
#define FCCGEOMETRYSERVICES_DD4HEPGEOMETRYSVC_H 1

// ATS
#include "DD4hepGeometryInterfaces/IDD4hepGeometrySvc.h"

//Gaudi
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/MsgStream.h"
// Core module
#include "CoreInterfaces/ServiceBase.h"

//DD4Hep
#include "DD4hep/LCDD.h"



/** @class DD4hepGeometrySvc DD4hepGeometrySvc.h FCCGeometryServices/DD4hepGeometrySvc.h
 
   This service provides the DD4hepGeometry to the ACTS package.
 
   @author julia.hrdinka@cern.ch
 */

namespace DD4hep {
    namespace Geometry {
        class DetElement;
    }
}

namespace Acts {

    class DD4hepGeometrySvc: public ServiceBase, virtual public IDD4hepGeometrySvc {

    public:

        /** default constructor */
        DD4hepGeometrySvc(const std::string& name, ISvcLocator* svc);
        /** destructor */
        virtual ~DD4hepGeometrySvc();
        /** retrieve interface ID */
        static const InterfaceID& interfaceID() { return IID_IDD4hepGeometrySvc; }
        /** initialize function */
        virtual StatusCode initialize();
        /** finalize function */
        virtual StatusCode finalize();
        /** Query the interfaces.
         /   Input: riid, Requested interface ID
         /          ppvInterface, Pointer to requested interface **/
        StatusCode queryInterface( const InterfaceID& riid, void** ppvInterface );
        /** receive DD4hep Geometry */
        virtual const DD4hep::Geometry::DetElement& worldDetElement() const;

    private:
    
        
        DD4hep::Geometry::DetElement                m_worldDetElement; //!<the world detector element containing the whole detector in DD4hep
        std::vector<std::string>                    m_xmlFileNames; //!<XML-file with the detector description
        MsgStream                                   m_log; //!<output stream
    };
}


inline const DD4hep::Geometry::DetElement& Acts::DD4hepGeometrySvc::worldDetElement() const { return m_worldDetElement; }


#endif //FCCGEOMETRYSERVICES_DD4HEPGEOMETRYSVC_H
