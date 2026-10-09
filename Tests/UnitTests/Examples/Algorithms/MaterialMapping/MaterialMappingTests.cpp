// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/BinnedSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"

using namespace Acts;
using namespace ActsExamples;

namespace {
class RecordingWriter final : public IMaterialWriter {
 public:
  void writeMaterial(const TrackingGeometryMaterial& material) override {
    ++calls;
    result = material;
  }
  unsigned int calls = 0;
  TrackingGeometryMaterial result;
};
class EmptyMaterialTracks final : public IAlgorithm {
 public:
  EmptyMaterialTracks() : IAlgorithm("EmptyMaterialTracks") {
    m_output.initialize("material_tracks");
  }

  ProcessCode execute(const AlgorithmContext& context) const override {
    m_output(context, {});
    return ProcessCode::SUCCESS;
  }

 private:
  WriteDataHandle<std::unordered_map<std::size_t, RecordedMaterialTrack>>
      m_output{this, "MaterialTracks"};
};
}  // namespace

BOOST_AUTO_TEST_CASE(FinalizedMaterialIsAvailableBeforeDestruction) {
  auto surface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30., 100.);
  surface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  surface->assignSurfaceMaterial(std::make_shared<ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D(
          {AxisSpec::DeferredEquidistant(1), AxisSpec::DeferredEquidistant(1)}),
      MappingType::Default, "cylinder"));
  BinnedSurfaceMaterialAccumulator::Config accumulatorConfig;
  accumulatorConfig.materialSurfaces = {surface.get()};
  MaterialMapper::Config mapperConfig;
  mapperConfig.assignmentFinder =
      std::make_shared<IntersectionMaterialAssigner>(
          IntersectionMaterialAssigner::Config{});
  mapperConfig.surfaceMaterialAccumulator =
      std::make_shared<BinnedSurfaceMaterialAccumulator>(accumulatorConfig);
  MaterialMapping::Config config;
  config.materialMapper = std::make_shared<MaterialMapper>(mapperConfig);
  auto writer = std::make_shared<RecordingWriter>();
  config.materialWriters = {writer};
  {
    auto mapping = std::make_shared<MaterialMapping>(config);
    BOOST_CHECK_THROW(mapping->material(), std::logic_error);
    Sequencer::Config sequenceConfig;
    sequenceConfig.events = 1;
    sequenceConfig.numThreads = 1;
    sequenceConfig.trackFpes = false;
    Sequencer sequencer(sequenceConfig);
    sequencer.addAlgorithm(std::make_shared<EmptyMaterialTracks>());
    sequencer.addAlgorithm(mapping);
    BOOST_CHECK_EQUAL(sequencer.run(), EXIT_SUCCESS);
    BOOST_CHECK_EQUAL(writer->calls, 1);
    BOOST_REQUIRE_EQUAL(mapping->material().keyedSurfaces.size(), 1);
    BOOST_CHECK(mapping->material().keyedSurfaces.at("cylinder").material ==
                writer->result.keyedSurfaces.at("cylinder").material);
  }
  BOOST_CHECK_EQUAL(writer->calls, 1);
  // An algorithm that never ran must not finalize or write during destruction.
  {
    MaterialMapping unused(config);
  }
  BOOST_CHECK_EQUAL(writer->calls, 1);
}
