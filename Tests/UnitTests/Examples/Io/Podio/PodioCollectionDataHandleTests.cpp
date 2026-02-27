// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <memory>

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/TrackCollection.h>
#include <podio/CollectionBase.h>

using namespace ActsExamples;

namespace ActsTests {

namespace {
static const auto logger =
    Acts::getDefaultLogger("PodioHandleTest", Acts::Logging::VERBOSE);
}  // namespace

BOOST_AUTO_TEST_SUITE(PodioCollectionDataHandleSuite)

BOOST_AUTO_TEST_CASE(EmulationCompatibility) {
  DummySequenceElement dummyElement;

  DataHandleBase::StateMapType state;
  WhiteBoard::AliasMapType aliases;

  PodioCollectionWriteHandle<edm4hep::TrackCollection> writeHandle{
      &dummyElement, "OutputTracks"};
  PodioCollectionReadHandle<edm4hep::TrackCollection> typedReadHandle{
      &dummyElement, "InputTracks"};
  ConsumeDataHandle<std::unique_ptr<podio::CollectionBase>> baseConsumeHandle{
      &dummyElement, "InputPodioWriter"};
  PodioCollectionReadHandle<edm4hep::MCParticleCollection> wrongTypedReadHandle{
      &dummyElement, "InputParticles"};

  writeHandle.initialize("tracks");
  typedReadHandle.initialize("tracks");
  baseConsumeHandle.initialize("tracks");
  wrongTypedReadHandle.initialize("tracks");

  writeHandle.emulate(state, aliases, *logger);

  BOOST_CHECK_NO_THROW(typedReadHandle.emulate(state, aliases, *logger));
  BOOST_CHECK_THROW(wrongTypedReadHandle.emulate(state, aliases, *logger),
                    SequenceConfigurationException);
  BOOST_CHECK_NO_THROW(baseConsumeHandle.emulate(state, aliases, *logger));
}

BOOST_AUTO_TEST_CASE(RuntimeTypedReadAndBaseConsume) {
  DummySequenceElement dummyElement;
  WhiteBoard wb;

  PodioCollectionWriteHandle<edm4hep::TrackCollection> writeHandle{
      &dummyElement, "OutputTracks"};
  PodioCollectionReadHandle<edm4hep::TrackCollection> typedReadHandle{
      &dummyElement, "InputTracks"};
  ConsumeDataHandle<std::unique_ptr<podio::CollectionBase>> baseConsumeHandle{
      &dummyElement, "InputPodioWriter"};

  writeHandle.initialize("tracks");
  typedReadHandle.initialize("tracks");
  baseConsumeHandle.initialize("tracks");

  edm4hep::TrackCollection tracks;
  writeHandle(wb, std::move(tracks));

  const auto& readTracks = typedReadHandle(wb);
  BOOST_CHECK_EQUAL(readTracks.size(), 0u);

  auto baseCollection = baseConsumeHandle(wb);
  BOOST_REQUIRE(baseCollection != nullptr);
  BOOST_CHECK(dynamic_cast<edm4hep::TrackCollection*>(baseCollection.get()) !=
              nullptr);
}

BOOST_AUTO_TEST_CASE(RuntimeTypeMismatchOnRead) {
  DummySequenceElement dummyElement;
  WhiteBoard wb;

  WriteDataHandle<std::unique_ptr<podio::CollectionBase>> baseWriteHandle{
      &dummyElement, "BaseWrite"};
  PodioCollectionReadHandle<edm4hep::TrackCollection> typedReadHandle{
      &dummyElement, "InputTracks"};

  baseWriteHandle.initialize("collections");
  typedReadHandle.initialize("collections");

  baseWriteHandle(wb, std::make_unique<edm4hep::MCParticleCollection>());

  BOOST_CHECK_THROW(typedReadHandle(wb), std::out_of_range);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
