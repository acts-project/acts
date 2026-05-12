#pragma once
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"

#include <memory>
#include <string>
#include <unordered_map>

#include <vecmem/memory/host_memory_resource.hpp>
#include <detray/detectors/odd_metadata.hpp>

namespace ActsExamples {

/// Converts Acts TrackingGeometry → detray detector in memory,
/// and writes the detray↔Acts geometry ID map to the whiteboard.
class ActsToDetrayDetectorAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The Acts tracking geometry to convert
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Name of the beampipe volume (detray requires it at index 0)
    std::string beampipeVolumeName = "Beampipe";
    /// Whiteboard key for the detray→Acts map
    std::string outputDetrayToActsMap = "detray-to-acts-map";
    /// Whiteboard key for the detray host detector
    std::string outputDetrayDetector  = "detray-detector";
    /// Optionally write detray JSON files
    std::string outputJsonDir = "";
  };

  explicit ActsToDetrayDetectorAlg(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode initialize() override;
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  mutable std::unordered_map<std::uint64_t, Acts::GeometryIdentifier> m_detrayToActsMap;

  WriteDataHandle<std::unordered_map<std::uint64_t, Acts::GeometryIdentifier>>
      m_outputDetrayToActsMap{this, "OutputDetrayToActsMap"};
};

}  // namespace ActsExamples