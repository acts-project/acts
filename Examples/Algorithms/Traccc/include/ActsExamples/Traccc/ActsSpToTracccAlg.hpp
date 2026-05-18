#pragma once
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <traccc/edm/spacepoint_collection.hpp>
#include <string>
#include <unordered_map>
#include <vecmem/memory/host_memory_resource.hpp>

namespace ActsExamples {

class ActsSpToTracccAlg final : public IAlgorithm {
 public:
  struct Config {
    std::string inputSpacePoints          = "spacepoints";
    std::string outputTracccSpacepoints   = "acts-traccc-spacepoints";
  };

  explicit ActsSpToTracccAlg(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  const Config& config() const { return m_cfg; }

  mutable vecmem::host_memory_resource m_mr;


 private:
  Config m_cfg;

  ReadDataHandle<SpacePointContainer>
      m_inputSpacePoints{this, "inputSpacePoints"};
  WriteDataHandle<traccc::edm::spacepoint_collection::host>
      m_outputTracccSpacepoints{this, "outputTracccSpacepoints"};
};

}  // namespace ActsExamples