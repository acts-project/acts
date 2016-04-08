///////////////////////////////////////////////////////////////////
// GeometryJsonWriter.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core module
#include "Algebra/StringConverters.h"
// Geometry module
#include "JsonWriters/GeometryJsonWriter.h"
#include "Detector/TrackingGeometry.h"
#include "Detector/TrackingVolume.h"
#include "Detector/Layer.h"
#include "Surfaces/Surface.h"
#include "Surfaces/RectangleBounds.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Volumes/CylinderVolumeBounds.h"
// output stream
#include <sstream>

DECLARE_TOOL_FACTORY(Acts::GeometryJsonWriter)

// constructor
Acts::GeometryJsonWriter::GeometryJsonWriter(const std::string& t, const std::string& n, const IInterface* p) : 
  Acts::AlgToolBase(t,n,p),
  m_outputFileName("TrackingGeometry.json"),
  m_outputPrecision(6),
  m_firstLayerWritten(false),
  m_layerCounter(0),
  m_firstVolumeWritten(false),
  m_suppressEmptyObjects(true),
  m_onlyWriteVolumes(false)
{
    
    declareInterface<Acts::IGeometryProcessor>(this);
    
    declareProperty("OutputFileName",  m_outputFileName);
    declareProperty("OutputPrecision", m_outputPrecision);
    declareProperty("SuppressEmptyObjects", m_suppressEmptyObjects);
    declareProperty("OnlyWriteVolumes", m_onlyWriteVolumes);
    
    m_extraVolumesToKeep = {"AtlasInnerSector"};
    m_nameTranslations = { {"AtlasInnerSector", "Inner Detector"}, {"InDet::Detectors::Pixel", "Pixel"}, 
                        {"InDet::Detectors::Pixel::NegativeEndcap","NegativeEndcap"}, {"InDet::Detectors::Pixel::Barrel","Barrel"}, 
                        {"InDet::Detectors::Pixel::PositiveEndcap","PositiveEndcap"}, {"InDet::Detectors::SCT::NegativeEndcap","NegativeEndcap"}, 
                        {"InDet::Detectors::SCT::Barrel","Barrel"}, {"InDet::Detectors::SCT::PositiveEndcap","PositiveEndcap"},
                        {"InDet::Containers::Container","SCT"}};
}

// destructor
Acts::GeometryJsonWriter::~GeometryJsonWriter()
{}


StatusCode Acts::GeometryJsonWriter::initialize()
{
    //Tool needs to be initialized
    if (!AlgToolBase::initialize()) return StatusCode::FAILURE;
    // open the file for writing
    m_outputFile.open(m_outputFileName.c_str());
    m_outputFile << std::setiosflags(std::ios::fixed);
    m_outputFile << std::setprecision(3);        
    // return the base::initialize() state
    return StatusCode::SUCCESS;    
}

StatusCode Acts::GeometryJsonWriter::finalize()
{
    // close the file
    m_outputFile.close();
    // return the base::finalize() state    
    return StatusCode::SUCCESS;    
}


StatusCode Acts::GeometryJsonWriter::process(const Acts::TrackingGeometry& tgeo) {
  
  MSG_VERBOSE("Start processing the TrackingGeometry recursively");
  // retrieve the highest tracking volume
  const Acts::TrackingVolume* worldVolume = tgeo.highestTrackingVolume();  
  if (worldVolume){
      MSG_VERBOSE("TrackingVolume '" << worldVolume->volumeName() << "' retrieved as highest level node.");
      return process(*worldVolume, 0);
  }
  // abort job
  MSG_FATAL("No highest level TrackingVolume found. Stopping recursive parsing, abort job.");
  return StatusCode::FAILURE;
}

// Processor Action to work on TrackingVolumes
StatusCode Acts::GeometryJsonWriter::process(const Acts::TrackingVolume& tvol, size_t level) {
  std::stringstream displayBuffer;
  for (size_t il = 0; il < level; ++il) displayBuffer << " ";
  MSG_VERBOSE(displayBuffer.str() << "TrackingVolume '" << tvol.volumeName() << "'");
  
  std::stringstream vstream;
  bool volumeHasMeaningfulContent;
  if (processNodeUncertain(tvol, vstream, volumeHasMeaningfulContent, ++level).isFailure() ){
    MSG_FATAL("Failed to call process(const TrackingVolume&) on confined volumes. Aborting.");
    return StatusCode::FAILURE;
  }
  MSG_VERBOSE("process Volume: Adding the following vstream to output file: "<<vstream.str());
  m_outputFile<<vstream.str();
  return StatusCode::SUCCESS; 
}

// Processor Action to work on Layers 
StatusCode Acts::GeometryJsonWriter::process(const Acts::Layer&, size_t) 
{
    // return SUCCESS
    return StatusCode::SUCCESS;
}

// Processor Action to work on Layers 
StatusCode Acts::GeometryJsonWriter::process(const Acts::Surface&, size_t) 
{
    // return SUCCESS
    return StatusCode::SUCCESS;
}

// Processor Action to work on TrackingVolumes
StatusCode Acts::GeometryJsonWriter::processNodeUncertain(const Acts::TrackingVolume& tvol, std::stringstream& stream, bool& hasMeaningfulContent,size_t level) {
  std::stringstream displayBuffer;
  for (size_t il = 0; il < level; ++il) displayBuffer << " ";
  MSG_VERBOSE(displayBuffer.str() << "processNodeUncertain: TrackingVolume '" << tvol.volumeName() << "'");


  hasMeaningfulContent = false;
  std::stringstream bstream; // info from the base volume object (i.e. this object)
  std::stringstream lstream; // info from layers
  std::stringstream vstream; // info from volumes
  bool firstcontent=true;

  // Process the contained layers if they exist
  const Acts::LayerArray* layerArray = tvol.confinedLayers();
  if (layerArray) {
    // display output
    auto& layers = layerArray->arrayObjects();
    MSG_VERBOSE(displayBuffer.str() << "--> has " << layers.size() << " confined layers." ); 

    bool layerHasMeaningfulContent = false;
    for (auto& layIter : layers){
      layerHasMeaningfulContent = false;
      if (!layIter){
        MSG_WARNING("Zero-pointer found in LayerArray - indicates problem !");
        continue;
      }
      std::stringstream templstream;
      if (processNodeUncertain(*layIter, templstream, layerHasMeaningfulContent, level).isFailure()){
        MSG_FATAL("Failed to call process(const Layer&) on confined layers. Aborting.");
        return StatusCode::FAILURE;
      }
      // MSG_VERBOSE(displayBuffer.str() <<"Processed layer and got back: "<<templstream.str());
      if (layerHasMeaningfulContent || !m_suppressEmptyObjects){
        if (firstcontent) {
          firstcontent=false;
        } else {
          lstream << ",";
        }
        lstream<<templstream.str();
        MSG_VERBOSE(displayBuffer.str() <<" -> Layer has meaningful content...");
        hasMeaningfulContent=true;
      }
    }
  } 

  if (!lstream.str().empty()) MSG_VERBOSE(displayBuffer.str() <<"Processed layers and lstream is "<<lstream.str());

  // Process the contained TrackingVolumes (recursively) if they exist
  auto confinedVolumes = tvol.confinedVolumes();
  // register the next round
  if (confinedVolumes) {
    bool volumeHasMeaningfulContent = false;
    auto& volumes = confinedVolumes->arrayObjects();
    MSG_VERBOSE(displayBuffer.str() << "--> has " << volumes.size() << " confined volumes." ); 
 
    for (auto& volumesIter : volumes){
      volumeHasMeaningfulContent = false;
      std::stringstream tempvstream;
     
      if (!volumesIter) {
        MSG_WARNING("Zero-pointer found in VolumeArray - indicates problem !");
        continue;
      }
     
      if (processNodeUncertain(*volumesIter, tempvstream, volumeHasMeaningfulContent, ++level).isFailure() ){
        MSG_FATAL("Failed to call process(const TrackingVolume&) on confined volumes. Aborting.");
        return StatusCode::FAILURE;
      }
     
      if (volumeHasMeaningfulContent || !m_suppressEmptyObjects) {
        if (firstcontent) {
          firstcontent=false;
        } else {
          vstream << ",";
        }
        vstream<<tempvstream.str();
        MSG_VERBOSE(displayBuffer.str() <<"Volume has meaningful content..."<<vstream.str());
        hasMeaningfulContent=true;
      }           
    }
  }

  if (!vstream.str().empty()) MSG_VERBOSE(displayBuffer.str() <<"Processed vols and vstream is "<<vstream.str());

  // Only write it out if it's interesting!
  if ( m_volumesToKeep.size()==0 || (m_volumesToKeep.size() && m_volumesToKeep.count(tvol.volumeName())) || m_extraVolumesToKeep.count(tvol.volumeName())){
    // So, volume is either in list, or the list is emtpy (and so we're keeping everything), or it's in the 'special' list.
    if (hasMeaningfulContent ) {
      if ( m_nameTranslations.count(tvol.volumeName()) ) {
        // Do we need to translate the name to something a bit nicer?
        stream << " { \"Name\" : \""<< m_nameTranslations[tvol.volumeName()]<<"\", ";
      } else {
        stream << " { \"Name\" : \""<< tvol.volumeName()<<"\", ";
      }
      dumpVolumeDetails(stream, tvol);
      stream <<", \"Layers\" : ["<< lstream.str()<<"],";
      stream <<" \"Volumes\" : ["<< vstream.str()<<"] ";
      stream << " }";
    } 
  } else {
    std::cout<<"m_volumesToKeep.size() "<<m_volumesToKeep.size()<<std::endl;
    MSG_VERBOSE(displayBuffer.str() <<"No interesting content or isn't on list of volumes to keep. ");
  }

  // return 
  return StatusCode::SUCCESS; 
}

void Acts::GeometryJsonWriter::dumpVolumeDetails(std::stringstream& stream, const Acts::TrackingVolume& tvol) const {
  dumpVolumeBounds(stream, tvol.volumeBounds());
  double cx = tvol.center().x();
  double cy = tvol.center().y();
  double cz = tvol.center().z();
  // get the euler angles         
  auto ea = tvol.transform().rotation().eulerAngles(0, 1, 2); 
  double e0 = ea[0];
  double e1 = ea[1]; 
  double e2 = ea[2];
  stream << ", \"Coords\": [ [" << cx << "," << cy << "," << cz << "],[" << e0 << "," << e1 << "," << e2 << "] ] ";
}

void Acts::GeometryJsonWriter::dumpVolumeBounds(std::stringstream& stream, const Acts::VolumeBounds& volumeBounds) const {
  const Acts::CylinderVolumeBounds* cvb = dynamic_cast<const Acts::CylinderVolumeBounds*>(&volumeBounds);
  if (cvb){
    stream << "\"Shape\" : \"CYL\",";
    stream << "\"Bounds\" : [" <<  cvb->innerRadius () << "," << cvb->outerRadius () <<"," << cvb->mediumRadius ()<<"," 
           << cvb->deltaRadius ()<<"," << cvb->halfPhiSector ()<<"," << cvb->halflengthZ()<<"]";
  }
  // If unknown, do nothing (for the moment).
}


StatusCode Acts::GeometryJsonWriter::processNodeUncertain(const Acts::Layer& lay, std::stringstream& stream, bool& hasMeangingfulContent, size_t level)
{    
    std::stringstream displayBuffer;
    for (size_t il = 0; il < level; ++il) displayBuffer << " ";
    
    MSG_VERBOSE(displayBuffer.str()<<"processNodeUncertainLayer: Dumping information for Layer with index " << lay.geoID().value());

    if (lay.surfaceArray()){
        size_t nSurfaces = lay.surfaceArray()->arrayObjects().size();
        // the layer has a surface Array - go for it
        // get the dimensions
        const Acts::Surface* referenceSurface = lay.surfaceArray()->arrayObjects()[0];
        // the reference Surface exists
        if (referenceSurface){
            // dynamic_cast to RectangleBounds
            const Acts::RectangleBounds* rBounds = dynamic_cast<const Acts::RectangleBounds*>(&(referenceSurface->bounds()));
            // we have rBounds - go on for the moment
            if (rBounds){
                hasMeangingfulContent = true;
                // the layer is defined
                stream << "{ ";
                stream << "\"Index\" :" << (lay.geoID().value()+(m_layerCounter++)) << ","; //!< @TODO remove layer counter
                stream << "\"Shape\" : \"BOX\",";
                stream << "\"Dimensions\" : [" <<  2*rBounds->halflengthX() << "," << 2.*rBounds->halflengthY() << ", 1.0 ], ";
                
                // count the surfaces
                size_t is = 1;
                // now loop of the surfaces and dump their position
                stream << "\"Coords\": [";
                for (auto& sf : lay.surfaceArray()->arrayObjects()){
                    // get x,y,z 
                    double cx = sf->center().x();
                    double cy = sf->center().y();
                    double cz = sf->center().z();
                    // get the euler angles         
                    auto ea = sf->transform().rotation().eulerAngles(0, 1, 2); 
                    double e0 = ea[0];
                    double e1 = ea[1]; 
                    double e2 = ea[2];
                    stream << "[" << cx << "," << cy << "," << cz << "],[" << e0 << "," << e1 << "," << e2 << "]";
                    if (is < nSurfaces) stream << ",";
                    ++is;
                }
                stream << " ]   }";
            }
            const Acts::TrapezoidBounds* tBounds = dynamic_cast<const Acts::TrapezoidBounds*>(&(referenceSurface->bounds()));
            if (tBounds){
              hasMeangingfulContent = true;
              // the layer is defined
              stream << "{ ";
              stream << "\"Index\" :" << (lay.geoID().value()+(m_layerCounter++)) << ",";
              stream << "\"Shape\" : \"TRA\",";
              stream << "\"Dimensions\" : [" <<  tBounds->minHalflengthX() << "," << tBounds->maxHalflengthX() <<","<< tBounds->halflengthY()<<"],";
              
              // count the surfaces
              size_t is = 1;
              // now loop of the surfaces and dump their position
              stream << "\"Coords\": [";
              for (auto& sf : lay.surfaceArray()->arrayObjects()){
                  // get x,y,z 
                  double cx = sf->center().x();
                  double cy = sf->center().y();
                  double cz = sf->center().z();
                  // get the euler angles         
                  auto ea = sf->transform().rotation().eulerAngles(0, 1, 2); 
                  double e0 = ea[0];
                  double e1 = ea[1]; 
                  double e2 = ea[2];
                  stream << "[" << cx << "," << cy << "," << cz << "],[" << e0 << "," << e1 << "," << e2 << "]";
                  if (is < nSurfaces) stream << ",";
                  ++is;
            }
            stream << " ]   }";
        }
        }
    }
    
    return StatusCode::SUCCESS;    
}


