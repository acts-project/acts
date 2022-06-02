#include "ExaTrkXTriton.hpp"

#include <iostream>

namespace tc = triton::client;

ExaTrkXTriton::ExaTrkXTriton(
    const std::string& modelName,
    const std::string& url,
    const std::string& modelVersion,
    uint32_t client_timeout, bool verbose
){
    options_ = std::make_unique<tc::InferOptions>(modelName);
    options_->model_version_ = modelVersion;
    options_->client_timeout_ = client_timeout;

    inputs_.clear();

    // Create a InferenceServerGrpcClient instance to communicate with the
    // server using gRPC protocol.
    std::unique_ptr<triton::client::InferenceServerGrpcClient> tClient;
    FAIL_IF_ERR(
        tc::InferenceServerGrpcClient::Create(&tClient, url, verbose),
        "unable to create grpc client");
    m_Client_ = std::move(tClient);
}

void ExaTrkXTriton::ClearInput()
{
  inputs_.clear();
}