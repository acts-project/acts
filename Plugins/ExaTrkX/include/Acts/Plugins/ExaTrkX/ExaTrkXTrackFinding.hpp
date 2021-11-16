#pragma once

#include <string>
#include <vector>

#include <core/session/onnxruntime_cxx_api.h>

class ExaTrkXTrackFinding {
  public:
    struct Config {
      /// input model directory
      std::string inputMLModuleDir;

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
    ExaTrkXTrackFinding(Config config);
    ExaTrkXTrackFinding();

    virtual ~ExaTrkXTrackFinding() {
    }

    void getTracks(
      std::vector<float>& input_values,
      std::vector<uint32_t>& spacepointIDs,
      std::vector<std::vector<uint32_t> >& trackCandidates) const;


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

  private:
    // configuration
    Config m_cfg;
    std::unique_ptr<Ort::Env> m_env;
    std::unique_ptr<Ort::Session> e_sess;
    std::unique_ptr<Ort::Session> f_sess;
    std::unique_ptr<Ort::Session> g_sess;
};