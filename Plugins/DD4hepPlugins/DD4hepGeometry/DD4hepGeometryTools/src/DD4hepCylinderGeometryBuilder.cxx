#include "DD4hepGeometryTools/DD4hepCylinderGeometryBuilder.h"
// Geometry module
#include "Detector/TrackingGeometry.h"
#include "Detector/TrackingVolume.h"
// DD4hepPlugin
#include "DD4hepGeometryUtils/DD4hepGeometryHelper.h"
// Core module
#include "CoreInterfaces/MsgBase.h"

DECLARE_COMPONENT(Add4hep::DD4hepCylinderGeometryBuilder)

Add4hep::DD4hepCylinderGeometryBuilder::DD4hepCylinderGeometryBuilder(const std::string& t,const std::string& n,const IInterface* p) :
Ats::AlgToolBase(t,n,p),
m_DD4hepGeometrySvc("", name),
m_volumeBuilder(nullptr),
m_DD4hepLayerHelper(nullptr),
m_trackingVolumeHelper(nullptr)
{
    declareInterface<ITrackingGeometryBuilder>(this);
    declareProperty("DD4hepGeometrySvc", m_DD4hepGeometrySvc);
    declareProperty("CylinderVolumeBuilder", m_volumeBuilder);
    declareProperty("DD4hepLayerHelper", m_DD4hepLayerHelper);
    declareProperty("CylinderVolumeHelper", m_trackingVolumeHelper);
}


Add4hep::DD4hepCylinderGeometryBuilder::~DD4hepCylinderGeometryBuilder()
{}

StatusCode Add4hep::DD4hepCylinderGeometryBuilder::initialize()
{
    // screen output in debug mode
    MSG_DEBUG( "initialize()" );
    //retrieve the Service providing the DD4hep geometry
    RETRIEVE_FATAL(m_DD4hepGeometrySvc);
    return StatusCode::SUCCESS;
}

StatusCode Add4hep::DD4hepCylinderGeometryBuilder::finalize()
{
    // screen output in debug mode
    MSG_DEBUG( "finalize()" );
    return StatusCode::SUCCESS;
}

Ats::TrackingGeometry* Add4hep::DD4hepCylinderGeometryBuilder::trackingGeometry() const
{
    // the return geometry -- and the highest volume
    Ats::TrackingGeometry* trackingGeometry = nullptr;
    Ats::TrackingVolumePtr    highestVolume = nullptr;
    Ats::TrackingVolumePtr beamPipeVolume   = nullptr;
    
    // get the DD4hep world detector element
    DD4hep::Geometry::DetElement detWorld = m_DD4hepGeometrySvc->worldDetElement();
    // get the sub detectors
    std::vector<DD4hep::Geometry::DetElement> detElements;
    const DD4hep::Geometry::DetElement::Children& children = detWorld.children();
    for (auto& detElement : children) detElements.push_back(detElement.second);
    //sort by id to build detector from bottom to top - new convention
    sort(detElements.begin(),detElements.end(),
         [](const DD4hep::Geometry::DetElement& a,
            const DD4hep::Geometry::DetElement& b) {
             return (a.id()<b.id());}
         );
    // loop over the volumes
    for (auto& detElement : detElements) {
        if (detElement.type()=="beamtube") {
            MSG_DEBUG("BeamPipe is being built");
            //extract material
            DD4hep::Geometry::Material mat = detElement.volume().material();
            //create the tracking volume
            beamPipeVolume = Ats::TrackingVolume::create(Add4hep::DD4hepGeometryHelper::extractTransform(detElement),Add4hep::DD4hepGeometryHelper::extractVolumeBounds(detElement),Ats::Material(mat.radLength(),mat.intLength(),mat.A(),mat.Z(),mat.density()),nullptr,nullptr,nullptr,nullptr,"BeamTube");
        }
        else {
        // assign a new highest volume (and potentially wrap around the given highest volume so far)
         highestVolume = m_volumeBuilder->trackingVolume(highestVolume,nullptr,m_DD4hepLayerHelper->createLayerTriple(detElement));
        }
    }
    // if you have a highest volume, stuff it into a TrackingGeometry
    if (highestVolume) {
        // see if the beampipe needs to be wrapped
        if (beamPipeVolume) highestVolume = m_trackingVolumeHelper->createContainerTrackingVolume({beamPipeVolume,highestVolume});
        
        // create the TrackingGeometry
        trackingGeometry = new Ats::TrackingGeometry(highestVolume);
    }
    // return the geometry to the service
    return trackingGeometry;
}