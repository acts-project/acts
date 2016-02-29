///////////////////////////////////////////////////////////////////
// GeometryJsonWriter.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef ATS_JSONWRITERS_GEOMETRYJSONDUMPER_H
#define ATS_JSONWRITERS_GEOMETRYJSONDUMPER_H

// Core module
#include "CoreInterfaces/AlgToolBase.h"
// Geometry module
#include "GeometryInterfaces/IGeometryProcessor.h"
// STL
#include <algorithm>
#include <string>
#include <fstream>


namespace Ats {

    class TrackingVolume;
    class Layer;
    class Surface;
    class VolumeBounds;

    /** @class GeometryJsonWriter

        The GeometryJsonWriter writes geometrical
        information of the tracking geomtetry 
        into an Ascii file for comparison

        @author Edward.Moyse@cern.ch, Andreas.Salzburger@cern.ch   
     */

    class GeometryJsonWriter : public AlgToolBase, virtual public IGeometryProcessor {

      public:

        /** Constructor */
        GeometryJsonWriter(const std::string&,const std::string&,const IInterface*);

        /** Destructor */
        virtual ~GeometryJsonWriter();

        /** AlgTool initialize method */
        StatusCode initialize();
        
        /** AlgTool finalize method */
        StatusCode finalize();

        /** Processor Action to work on TrackingGeometry& tgeo */
        virtual StatusCode process(const TrackingGeometry& tgeo);
       
        /** Processor Action to work on TrackingVolumes - the level is for the hierachy tree*/
        virtual StatusCode process(const TrackingVolume& tvol, size_t level = 0);   
       
        /** Processor Action to work on Layers */
        virtual StatusCode process(const Layer& lay, size_t level = 0);
       
        /** Processor Action to work on Surfaces */
        virtual StatusCode process(const Surface& surf, size_t level = 0);

      private:

        /** Current implementation: write root visualization to file stream */
        StatusCode processNodeUncertain(const Layer& lay, std::stringstream& stream, bool& hasMeangingfulContent, size_t level);
        StatusCode processNodeUncertain(const TrackingVolume& tvol, std::stringstream& stream, bool& hasMeangingfulContent, size_t level);

        /** Current implementation: write root visualization to file stream */
        void dumpVolumeDetails(std::stringstream& stream, const TrackingVolume& tvol) const ;
        void dumpVolumeBounds(std::stringstream& stream, const VolumeBounds& volumeBounds) const ;

        std::ofstream                       m_outputFile;
        std::string                         m_outputFileName;  //!< where the tree is written to  
        int                                 m_outputPrecision;
        bool                                m_firstLayerWritten;
        int                                 m_layerCounter;
        bool                                m_firstVolumeWritten;
        bool                                m_suppressEmptyObjects; //!< If true, don't write out empty layers and volumes.
        bool                                m_onlyWriteVolumes; //!< If true, only write out the volumes.
        std::map<std::string, std::string>  m_nameTranslations; //!< If we want to have simplified names for layers, volumes etc. Set in initialise().
        std::set<std::string>               m_layersToKeep; //!< if empty, we take everything. If not, we only take layers listed here. Set in initialise().
        std::set<std::string>               m_volumesToKeep; //!< if empty, we take everything that is 'interesting'. If not, we will always take volumes listed here. Set in initialise().
        std::set<std::string>               m_extraVolumesToKeep; //!< Any volumes named here will always be kept, whether they pass other criteria or not. Set in initialise().
 
    };
}

#endif

