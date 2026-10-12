// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include "ReconstructionFixture.hpp"

namespace ActsTests {
struct CalorimeterIoEvent {
  ActsExamples::WhiteBoard board;
  ActsExamples::AlgorithmContext context{0, 0, board, 0};
  DummySequenceElement source;
  explicit CalorimeterIoEvent(
      podio::Frame frame = ActsExamples::CalorimeterFixture::makeFrame(true)) {
    write("events", std::move(frame));
  }
  template <typename T>
  void write(const std::string& name, T data) {
    ActsExamples::WriteDataHandle<T> handle(&source, "Write");
    handle.initialize(name);
    handle(context, std::move(data));
  }
  template <typename T>
  const T& read(const std::string& name) {
    ActsExamples::ReadDataHandle<T> handle(&source, "Read");
    handle.initialize(name);
    return handle(context);
  }
  template <typename T>
  void writePodio(const std::string& name, T data) {
    ActsExamples::PodioCollectionWriteHandle<T> handle(&source, "WritePodio");
    handle.initialize(name);
    handle(context, std::move(data));
  }
  template <typename T>
  const T& readPodio(const std::string& name) {
    ActsExamples::PodioCollectionReadHandle<T> handle(&source, "ReadPodio");
    handle.initialize(name);
    return handle(context);
  }
  void reconstruct(double energyScale = 2) {
    namespace Fixture = ActsExamples::CalorimeterFixture;
    ActsExamples::EDM4hepCalorimeterInputConverter{Fixture::inputConfig()}
        .execute(context);
    auto response = Fixture::responseConfig();
    response.response.energyScale = energyScale;
    ActsExamples::CalorimeterDigitizationAlgorithm{response}.execute(context);
    ActsExamples::CalorimeterClusteringAlgorithm{Fixture::clusteringConfig()}
        .execute(context);
    ActsExamples::EDM4hepCalorimeterOutputConverter{Fixture::hitOutputConfig()}
        .execute(context);
  }
  void persistClusters() {
    ActsExamples::EDM4hepCalorimeterClusterOutputConverter{
        ActsExamples::CalorimeterFixture::clusterOutputConfig()}
        .execute(context);
  }
};
inline podio::Frame emptyCalorimeterFrame() {
  podio::Frame frame;
  frame.put(edm4hep::SimCalorimeterHitCollection{}, "SimCaloHits");
  return frame;
}
}  // namespace ActsTests
