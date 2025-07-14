#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <string>
#include <limits>
#include <vector>

namespace ActsExamples {

class MeasurementCutter final : public IAlgorithm {
    public:
    struct Config {
        std::string inputSpacePoints;
        std::string inputMeasurementParticlesMap;
        std::string inputParticles;
        
        std::string outputParticleMeasurementsMap;
        std::string outputMeasurementsParticlesMap;
        std::string outputParticles;

        // Support both single values and arrays for R bounds
        std::vector<double> maxR = {std::numeric_limits<double>::infinity()};
        std::vector<double> maxZ = {std::numeric_limits<double>::infinity()};
        
        double minR = 0;
        double minZ = -std::numeric_limits<double>::infinity();
    };

    MeasurementCutter(const Config& config, Acts::Logging::Level level);

    ProcessCode execute(const AlgorithmContext& ctx) const final; 

    const Config& config() const { return m_cfg; }

    private:
        Config m_cfg;
        
        // Helper method to check if a point passes the selection criteria
        bool passesSelection(const SimSpacePoint& spacepoint) const;

    ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this, "input_space_points"};
    ReadDataHandle<IndexMultimap<ActsFatras::Barcode>> m_inputMeasurementParticlesMap{this, "input_measurement_particles_map"};
    ReadDataHandle<SimParticleContainer> m_inputParticles{this, "input_particles"};

    WriteDataHandle<InverseMultimap<ActsFatras::Barcode>> m_outputParticleMeasurementsMap{this, "output_particle_measurements_map"};
    WriteDataHandle<IndexMultimap<ActsFatras::Barcode>> m_outputMeasurementsParticlesMap{this, "output_measurements_particles_map"};
    WriteDataHandle<SimParticleContainer> m_outputParticles{this, "output_particles"};
};

}  // namespace ActsExamples