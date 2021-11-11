#pragma once

#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Index.hpp"

#include <string>
#include <vector>

#include <core/session/onnxruntime_cxx_api.h>

namespace ActsExamples {

class TrackFindingMLBasedAlgorithm final : public BareAlgorithm {
  public:
    struct Config {
      /// Input spacepoints collection.
      std::string inputSpacePoints;

      /// input model directory
      std::string inputMLModuleDir;

      /// Output protoTracks collection.
      std::string outputProtoTracks;

      // hyperparameters in the pipeline.
      int64_t spacepointFeatures = 3;
      int embeddingDim = 8;
      float rVal = 1.6;
      int knnVal = 500;
      float filterCut = 0.21;
    };

    /// Constructor of the track finding algorithm
    ///
    /// @param cfg is the config struct to configure the algorithm
    /// @param level is the logging level
    TrackFindingMLBasedAlgorithm(Config cfg, Acts::Logging::Level lvl);

    virtual ~TrackFindingMLBasedAlgorithm() {
      delete e_sess;
      delete f_sess;
      delete g_sess;
      delete m_env;
    }

    /// Framework execute method of the track finding algorithm
    ///
    /// @param ctx is the algorithm context that holds event-wise information
    /// @return a process code to steer the algorithm flow
    ActsExamples::ProcessCode execute(
        const ActsExamples::AlgorithmContext& ctx) const final;

    const Config& config() const { return m_cfg; }

  private:
    void initTrainedModels();
    
    void runSessionWithIoBinding(
      Ort::Session& sess,
      std::vector<const char*>& inputNames,
      std::vector<Ort::Value> & inputData,
      std::vector<const char*>& outputNames,
      std::vector<Ort::Value>&  outputData) const;

    void buildEdges(std::vector<float>& embedFeatures,
          std::vector<int64_t>& edgeList, int64_t numSpacepoints) const;

    void getTracks(
      std::vector<float>& input_values,
      std::vector<ActsExamples::Index>& spacepointIDs,
      ProtoTrackContainer& trackCandidates) const;

  private:
    // configuration
    Config m_cfg;
    Ort::Env* m_env;
    Ort::Session* e_sess;
    Ort::Session* f_sess;
    Ort::Session* g_sess;

};

}