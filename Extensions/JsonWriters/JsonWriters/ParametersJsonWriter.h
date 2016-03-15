///////////////////////////////////////////////////////////////////
// GeometryJsonWriter.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_JSONWRITERS_PARAMETERSBASEJSONEVENTWRITER_H
#define ACTS_JSONWRITERS_PARAMETERSBASEJSONEVENTWRITER_H 1

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Test module
#include "TestInterfaces/IParametersBaseProcessor.h"
// EventData module
#include "ParametersBase/TrackParametersBase.h"
// STL
#include <algorithm>
#include <string>
#include <fstream>


namespace Acts {

    class TrackingVolume;
    class Layer;
    class Surface;
    class VolumeBounds;

    /** @class ParametersJsonWriter

        write track parameters into a Json file

        @author Edward.Moyse@cern.ch, Andreas.Salzburger@cern.ch   
     */

    class ParametersJsonWriter : public AlgToolBase, virtual public IParametersBaseProcessor {

      public:

        /** Constructor */
        ParametersJsonWriter(const std::string&,const std::string&,const IInterface*);

        /** Destructor */
        virtual ~ParametersJsonWriter();

        /** AlgTool initialize method */
        StatusCode initialize();
        
        /** AlgTool finalize method */
        StatusCode finalize();

        /** Processor Action to work on ParameterBase   */
        StatusCode process(const TrackParametersBase& parameters);
        
        /** Processor Action to work on ParameterBase vector  */
        StatusCode process(const std::vector<const TrackParametersBase*>& pBaseVector);

        /** Processor Action initializations */
        StatusCode initProcessor();

      private:
        StatusCode openFile(); 
        StatusCode closeFile(); 

        std::ofstream         m_outputFile;
        std::string           m_outputFileBase;  
        std::string           m_outputType;
        std::string           m_outputName;
        std::string           m_outputObject;
        int                   m_outputPrecision;
        int                   m_runNumber;
        int                   m_eventNumber;
        int                   m_trackNumber;
       
    };
}

#endif // ACTS_JSONWRITERS_PARAMETERSBASEJSONEVENTWRITER_H

