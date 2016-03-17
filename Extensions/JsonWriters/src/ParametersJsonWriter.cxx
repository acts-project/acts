///////////////////////////////////////////////////////////////////
// ParametersJsonWriter.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core module
#include "Algebra/StringConverters.h"
// Geometry module
#include "JsonWriters/ParametersJsonWriter.h"
// Event Data Module
#include "ParametersBase/TrackParametersBase.h"
// output stream
#include <sstream>
#include <boost/lexical_cast.hpp>

DECLARE_TOOL_FACTORY(Acts::ParametersJsonWriter)

// constructor
Acts::ParametersJsonWriter::ParametersJsonWriter(const std::string& t, const std::string& n, const IInterface* p) : 
  Acts::AlgToolBase(t,n,p),
  m_outputFileBase("EventDump_"),
  m_outputType("xAOD::Type::TrackParticle"),
  m_outputName("InDetTrackParticles"),
  m_outputObject("Trk "),
  m_outputPrecision(3),
  m_runNumber(0),
  m_eventNumber(0),
  m_trackNumber(0)

{
    declareInterface<Acts::IParametersBaseProcessor>(this);
    
    declareProperty("OutputFileBase",  m_outputFileBase);
    declareProperty("OutputType",      m_outputType);
    declareProperty("OutputName",      m_outputName);
    declareProperty("OutputObject",    m_outputObject);
    declareProperty("OutputPrecision", m_outputPrecision);
    
    declareProperty("RunNumber",       m_runNumber);
    declareProperty("EventNumber",     m_eventNumber);
}  

// destructor
Acts::ParametersJsonWriter::~ParametersJsonWriter()
{}


StatusCode Acts::ParametersJsonWriter::initialize()
{
    //Tool needs to be initialized
    if (!AlgToolBase::initialize()) return StatusCode::FAILURE;
    return openFile();    
}

StatusCode Acts::ParametersJsonWriter::finalize()
{
    return closeFile();
}

StatusCode Acts::ParametersJsonWriter::process(const TrackParametersBase& tpBase) 
{

    std::string trackObject = m_outputObject;
    trackObject += boost::lexical_cast<std::string>(m_trackNumber);
    
    if (m_trackNumber) m_outputFile << ","; // from the previous track
    m_outputFile << "\"" << trackObject << "\": {";
    m_outputFile << "\"chi2\": 12.6197, \"dof\": 15,";
    m_outputFile << "\"dparams\": ["  << tpBase.parameters()[eLOC_D0] 
                              << ", " << tpBase.parameters()[eLOC_Z0]
                              << " ," << tpBase.parameters()[ePHI] 
                              << ", " << tpBase.parameters()[eTHETA] 
                              << " ," << tpBase.parameters()[eQOP] << " ],";
    
    ++m_trackNumber;
    
    return StatusCode::SUCCESS;    
}

StatusCode Acts::ParametersJsonWriter::process(const std::vector<const TrackParametersBase*>& pBaseVector) {
    
    // prepare
    m_outputFile << "\"pos\": [";
    for (size_t ipb = 0; ipb < pBaseVector.size(); ++ipb){
        // get the position
        const Vector3D& pos = pBaseVector[ipb]->position();
        m_outputFile << "[" << pos.x() << ", " << pos.y() << ", " << pos.z() << "]";
        if ( (ipb+1) != pBaseVector.size() ) m_outputFile << ",";
    }  
    m_outputFile << "]"; // closes the pos array
    m_outputFile << "}"; // closes teh Output object
    
    return StatusCode::SUCCESS;
    
}

StatusCode Acts::ParametersJsonWriter::initProcessor() 
{
    if (closeFile().isFailure())
        MSG_WARNING("Problem with file opening.");
    if (openFile().isFailure())
        MSG_WARNING("Problem with file opening.");
    m_trackNumber = 0;
    ++m_eventNumber;
    
    return StatusCode::SUCCESS;
}


StatusCode Acts::ParametersJsonWriter::openFile() {
    // open the file for writing
    std::string outputFileNameEvent = m_outputFileBase;
    outputFileNameEvent += boost::lexical_cast<std::string>(m_eventNumber);
    outputFileNameEvent += ".json";
    // open a new file 
    m_outputFile.open(outputFileNameEvent.c_str());
    //m_outputFile << std::setiosflags(std::ios::fixed);
    //m_outputFile << std::setprecision(m_outputPrecision); 
    // 
    m_outputFile << "{ \"event number\":" << m_eventNumber << ", \"run number\": 0, ";
    m_outputFile << "\"" << m_outputType << "\" : {"; 
    m_outputFile << "\"" << m_outputName << "\" : {";
    
    return StatusCode::SUCCESS;    
}

StatusCode Acts::ParametersJsonWriter::closeFile() 
{
    
    // close the json format
    m_outputFile << "}";  // closes outptu name 
    m_outputFile << "}";  // closes output type
    m_outputFile << "}";  // closes event
    // close the file
    m_outputFile.close();
    // return the base::finalize() state    
    return StatusCode::SUCCESS;        
}

