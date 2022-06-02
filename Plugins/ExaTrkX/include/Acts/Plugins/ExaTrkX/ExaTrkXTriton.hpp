#pragma once
#include <torch/torch.h>
#include <torch/script.h>
using namespace torch::indexing;


#include "grpc_client.h"
#include "grpc_service.pb.h"

#include <string>
#include <vector>
#include <memory>

namespace tc = triton::client;


#define FAIL_IF_ERR(X, MSG)                                        \
{                                                                \
    tc::Error err = (X);                                          \
    if (!err.IsOk()) {                                             \
      std::cerr << "error: " << (MSG) << ": " << err << std::endl; \
      exit(1);                                                     \
    }                                                              \
}

class ExaTrkXTriton {
  public:
    ExaTrkXTriton(const std::string& modelName,
    const std::string& url,
    const std::string& modelVersion="",
    uint32_t client_timeout=0, bool verbose=false);


    ExaTrkXTriton() = delete;
    ExaTrkXTriton(const ExaTrkXTriton&) = delete;
    ExaTrkXTriton& operator=(const ExaTrkXTriton&) = delete;
    ~ExaTrkXTriton() {};

    template <typename T>
    bool AddInput(
      const std::string& inputName, const std::vector<int64_t>& inputShape,
      std::vector<T>& inputValues){
      tc::InferInput* input0;
      std::string dataType{"FP32"};
      if (std::is_same<T, int32_t>::value) {
        dataType = "INT32";
      } else if (std::is_same<T, int64_t>::value) {
        dataType = "INT64";
      } else {}

      FAIL_IF_ERR(
          tc::InferInput::Create(&input0, inputName, inputShape, dataType), 
          "unable to get"+inputName);

      std::shared_ptr<tc::InferInput> input0_ptr(input0);

      FAIL_IF_ERR(input0_ptr->AppendRaw(
          reinterpret_cast<uint8_t*>(&inputValues[0]), // why uint8?
          inputValues.size() * sizeof(T)), "unable to set data"+inputName);

      inputs_.push_back(input0_ptr);

      return true;
    }

    template <typename T>
    bool AddInputFromTorch(const std::string& inputName, const at::Tensor& inputTensor){
      std::vector<int64_t> inputShape{inputTensor.sizes().vec()};  
      std::vector<T> inputValues;
      std::copy(
        inputTensor.data_ptr<T>(),
        inputTensor.data_ptr<T>() + inputTensor.numel(),
        std::back_inserter(inputValues));
      
      return AddInput<T>(inputName, inputShape, inputValues);
    }

    void ClearInput();

    template <typename T>
    bool GetOutput(const std::string& outputName,
      std::vector<T>& outputData, const std::vector<int64_t>&  outputShape)
    {
      // std::cout << "In the inference" << std::endl;
      if (outputShape.size() != 2) {
          std::cerr << "error: output shape must be 2D" << std::endl;
      }
      tc::Headers http_headers;
      grpc_compression_algorithm compression_algorithm =
          grpc_compression_algorithm::GRPC_COMPRESS_NONE;

      tc::InferResult* results;
      std::vector<const tc::InferRequestedOutput*> outputs = {};

      // std::cout << "prepare to run inference with " << inputs_.size() << " input(s)." << std::endl;
      std::vector<tc::InferInput*> inputs;
      for (auto& input : inputs_) {
          inputs.push_back(input.get());
      }
      FAIL_IF_ERR(m_Client_->Infer(
          &results, *options_, inputs, outputs, http_headers,
          compression_algorithm), "unable to run inference");

      std::shared_ptr<tc::InferResult> results_ptr;
      results_ptr.reset(results);

      T* output_data;
      size_t output_size;
      results_ptr->RawData(
            outputName, (const uint8_t**)&output_data, &output_size);

      // std::cout << "before transfer data: " << output_size << " " << sizeof(output_data) << std::endl;
      outputData.clear();
      int64_t output_entries = output_size / sizeof(output_data);
      // std::cout << "output_size: " << output_entries << " " << output_data[0] << " " << output_data[1] << std::endl;

      int64_t num_rows = outputShape[0];
      int64_t num_cols = outputShape[1];
      if (num_rows < 0) {
          num_rows = (int64_t) output_entries / num_cols;
      }
      if (num_cols < 0) {
          num_cols = (int64_t) output_entries / num_rows;
      }

      // std::cout << "output_shape: " << num_rows << " " << num_cols << std::endl;
      for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) 
          outputData.push_back(output_data[i*num_cols + j]);
      }
      return true;
    }
  
  private:
    std::unique_ptr<tc::InferenceServerGrpcClient> m_Client_;
    std::vector<std::shared_ptr<tc::InferInput> > inputs_;
    std::unique_ptr<tc::InferOptions> options_;
};
