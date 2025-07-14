#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"




#include "ActsExamples/TruthTracking/MeasurementCutter.hpp"



namespace ActsExamples {


MeasurementCutter::MeasurementCutter(const Config& config, Acts::Logging::Level level) 
    : IAlgorithm("MeasurementCutter", level), m_cfg(config) {

        // Validate array bounds configuration
        if (m_cfg.maxR.size() != m_cfg.maxZ.size()) {
            throw std::invalid_argument("maxR and maxZ arrays must have the same length");
        }
        
        // Validate single value bounds
        if (m_cfg.minR < 0) {
            throw std::invalid_argument("minR must be non-negative");
        }
        
        // Validate input/output configuration
        if (m_cfg.inputSpacePoints.empty()) {
            throw std::invalid_argument("Missing input space points");
        }
        if (m_cfg.inputMeasurementParticlesMap.empty()) {
            throw std::invalid_argument("Missing input measurement particles map");
        }
        if (m_cfg.outputParticleMeasurementsMap.empty()) {
            throw std::invalid_argument("Missing output particle measurements map");
        }
        if (m_cfg.outputMeasurementsParticlesMap.empty()) {
            throw std::invalid_argument("Missing output measurements particles map");
        }
        if (m_cfg.inputParticles.empty()) {
            throw std::invalid_argument("Missing input particles");
        }
        if (m_cfg.outputParticles.empty()) {
            throw std::invalid_argument("Missing output particles");
        }

        m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
        m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
        m_inputParticles.initialize(m_cfg.inputParticles);

        m_outputParticleMeasurementsMap.initialize(m_cfg.outputParticleMeasurementsMap);
        m_outputMeasurementsParticlesMap.initialize(m_cfg.outputMeasurementsParticlesMap);
        m_outputParticles.initialize(m_cfg.outputParticles);
    }

bool MeasurementCutter::passesSelection(const SimSpacePoint& spacepoint) const {
    const auto r = spacepoint.r();
    const auto z = spacepoint.z();
    
    // Apply minimum bounds
    if (r < m_cfg.minR || z < m_cfg.minZ) {
        return false;
    }
    
    // Apply array-based bounds (similar to graph_modifier logic)
    bool passesArrayBounds = false;
    for (size_t i = 0; i < m_cfg.maxR.size(); ++i) {
        if (r <= m_cfg.maxR[i] && z <= m_cfg.maxZ[i]) {
            passesArrayBounds = true;
            break;
        }
    }
    
    return passesArrayBounds;
}

ProcessCode MeasurementCutter::execute(const AlgorithmContext& ctx) const {
    const SimSpacePointContainer& spacePoints = m_inputSpacePoints(ctx);
    const IndexMultimap<ActsFatras::Barcode>& measurementParticlesMap = m_inputMeasurementParticlesMap(ctx);
    const SimParticleContainer& particles = m_inputParticles(ctx);

    ACTS_VERBOSE("Total space points available: " << spacePoints.size());
    ACTS_VERBOSE("Total measurements in map: " << measurementParticlesMap.size());

    SimParticleContainer outputParticles;
    outputParticles.reserve(particles.size());

    IndexMultimap<ActsFatras::Barcode> outputMeasurementsParticlesMap;
    outputMeasurementsParticlesMap.reserve(measurementParticlesMap.size());
    
    std::set<Index> processedMeasurements; // Track which measurements we've already processed
    
    size_t foundSpacePoints = 0;
    size_t missingSpacePoints = 0;
    size_t extendable = 0; // If the track can be extended by the CKF
    size_t singleMeasurementExtendable = 0; // If the track can be extended by a single measurement
    
    // First loop: Process spacepoints and their associated measurements
    for (const auto& spacepoint : spacePoints) {
        // Apply selection criteria using the helper method
        bool keepSpacePoint = passesSelection(spacepoint);
        
        if (keepSpacePoint) {
            ACTS_VERBOSE("Selected SpacePoint: " << spacepoint.x() << " " << spacepoint.y() << " " << spacepoint.z() << " " << spacepoint.r());
        } else {
            ACTS_DEBUG("Removed SpacePoint: " << spacepoint.x() << " " << spacepoint.y() << " " << spacepoint.z() << " " << spacepoint.r());
            extendable += spacepoint.sourceLinks().size();
        }
        
        
        // Process all source links for this spacepoint
        for (const auto& sl : spacepoint.sourceLinks()) {
            Index measIndex = sl.template get<IndexSourceLink>().index();
            ACTS_VERBOSE("Processing source link measurement index: " << measIndex);
            
            // Find the particle ID for this measurement
            auto range = measurementParticlesMap.equal_range(measIndex);
            for (auto it = range.first; it != range.second; ++it) {
                const ActsFatras::Barcode& particleId = it->second;
                foundSpacePoints++;
                
                //ACTS_DEBUG("Measurement " << measIndex << " with particle ID " << particleId
                //          << " has space point at (" << spacepoint.x() << ", "
                //          << spacepoint.y() << ", " << spacepoint.z() << ") r = "
                //          << spacepoint.r() << " | keep: " << std::boolalpha
                //          << keepSpacePoint << " | extendable: " << isExtendable << " | Source links: "
                //          << spacepoint.sourceLinks().size());
                
                if (keepSpacePoint) {
                    outputMeasurementsParticlesMap.insert({measIndex, particleId});
                }
                
                processedMeasurements.insert(measIndex);
            }
        }
        
        
    }
    
    // Second loop: Process measurements that don't have corresponding spacepoints
    for (const auto& [measIndex, particleId] : measurementParticlesMap) {
        if (processedMeasurements.find(measIndex) == processedMeasurements.end()) {
            missingSpacePoints++;
            // Measurement without space point - decide policy here
            ACTS_DEBUG("Measurement " << measIndex << " has no corresponding space point - removing");
            singleMeasurementExtendable++;
        }
    }
    
    ACTS_VERBOSE("Space points found for measurements: " << foundSpacePoints);
    ACTS_VERBOSE("Measurements without space points: " << missingSpacePoints);
    ACTS_DEBUG("Selected " << outputMeasurementsParticlesMap.size() << " out of " << measurementParticlesMap.size() << " measurements");
    
    // For now, copy all particles (you might want to add particle-level filtering later)
    for (const auto& particle : particles) {
        outputParticles.insert(particle);
    }

    ACTS_INFO("Extendable Spacepoint measurements: " << extendable);
    ACTS_INFO("Single strip measurement extendable measurements: " << singleMeasurementExtendable);

    m_outputParticleMeasurementsMap(ctx, invertIndexMultimap(outputMeasurementsParticlesMap));
    m_outputMeasurementsParticlesMap(ctx, std::move(outputMeasurementsParticlesMap));
    m_outputParticles(ctx, std::move(outputParticles));

    return ProcessCode::SUCCESS;
}


} // namespace ActsExamples

