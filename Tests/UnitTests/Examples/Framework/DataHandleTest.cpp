// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_log.hpp>

#include "Acts/Tests/CommonHelpers/WhiteBoardUtilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

using namespace Acts::Test;
using namespace ActsExamples;
using Acts::Logging::ScopedFailureThreshold;

// Global logger instance for all tests
const Acts::Logger& logger() {
  static const auto logger =
      Acts::getDefaultLogger("DataHandleTest", Acts::Logging::VERBOSE);
  return *logger;
}

BOOST_AUTO_TEST_SUITE(DataHandleTest)

BOOST_AUTO_TEST_CASE(BasicOperations) {
  WhiteBoard wb;
  DummySequenceElement dummyElement;

  BOOST_TEST_CHECKPOINT("Test adding and retrieving objects");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("test_key");
    writeHandle(wb, 42);

    ReadDataHandle<int> readHandle(&dummyElement, "test");
    readHandle.initialize("test_key");
    BOOST_CHECK_EQUAL(readHandle(wb), 42);
  }

  BOOST_TEST_CHECKPOINT("Test initialize with empty key");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    BOOST_CHECK_THROW(writeHandle.initialize(""), std::invalid_argument);
  }

  BOOST_TEST_CHECKPOINT("Test maybeInitialize");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.maybeInitialize("");  // Should not throw
    BOOST_CHECK(!writeHandle.isInitialized());

    writeHandle.maybeInitialize("maybe_key");  // Should initialize
    BOOST_CHECK(writeHandle.isInitialized());
    writeHandle(wb, 42);

    ReadDataHandle<int> readHandle(&dummyElement, "test");
    readHandle.maybeInitialize("maybe_key");
    BOOST_CHECK(readHandle.isInitialized());
    BOOST_CHECK_EQUAL(readHandle(wb), 42);
  }

  BOOST_TEST_CHECKPOINT("Test uninitialized handles");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    BOOST_CHECK_THROW(writeHandle(wb, 42), std::runtime_error);

    ReadDataHandle<int> readHandle(&dummyElement, "test");
    BOOST_CHECK_THROW(readHandle(wb), std::runtime_error);
  }

  BOOST_TEST_CHECKPOINT("Test adding duplicate objects");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("duplicate_key");
    writeHandle(wb, 42);
    BOOST_CHECK_THROW(writeHandle(wb, 43), std::invalid_argument);
  }

  BOOST_TEST_CHECKPOINT("Test getting non-existent objects");
  {
    ReadDataHandle<int> readHandle(&dummyElement, "test");
    readHandle.initialize("missing_key");
    BOOST_CHECK_THROW(readHandle(wb), std::out_of_range);
  }

  BOOST_TEST_CHECKPOINT("Test getting objects with wrong type");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("type_key");
    writeHandle(wb, 42);

    ReadDataHandle<std::string> readHandle(&dummyElement, "test");
    readHandle.initialize("type_key");
    BOOST_CHECK_THROW(readHandle(wb), std::out_of_range);
  }

  BOOST_TEST_CHECKPOINT("Test similar name suggestions");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("similar_key");
    writeHandle(wb, 42);

    ReadDataHandle<int> readHandle(&dummyElement, "test");
    readHandle.initialize("similr_key");  // Typo in key
    BOOST_CHECK_THROW(readHandle(wb), std::out_of_range);
  }
}

BOOST_AUTO_TEST_CASE(ComplexTypes) {
  WhiteBoard wb;
  DummySequenceElement dummyElement;

  BOOST_TEST_CHECKPOINT("Test with vector type");
  {
    std::vector<int> data = {1, 2, 3};
    WriteDataHandle<std::vector<int>> writeHandle(&dummyElement, "test");
    writeHandle.initialize("vector_key");
    writeHandle(wb, std::move(data));

    ReadDataHandle<std::vector<int>> readHandle(&dummyElement, "test");
    readHandle.initialize("vector_key");
    const auto& result = readHandle(wb);
    BOOST_CHECK_EQUAL(result.size(), 3);
    BOOST_CHECK_EQUAL(result[0], 1);
    BOOST_CHECK_EQUAL(result[1], 2);
    BOOST_CHECK_EQUAL(result[2], 3);
  }

  BOOST_TEST_CHECKPOINT("Test with string type");
  {
    std::string data = "test string";
    WriteDataHandle<std::string> writeHandle(&dummyElement, "test");
    writeHandle.initialize("string_key");
    writeHandle(wb, std::move(data));

    ReadDataHandle<std::string> readHandle(&dummyElement, "test");
    readHandle.initialize("string_key");
    BOOST_CHECK_EQUAL(readHandle(wb), "test string");
  }
}

BOOST_AUTO_TEST_CASE(DataHandleCompatibility) {
  WhiteBoard wb;
  DummySequenceElement dummyElement;

  BOOST_TEST_CHECKPOINT("Test write handle with same type");
  {
    WriteDataHandle<int> writeHandle1(&dummyElement, "test1");
    writeHandle1.initialize("same_key");
    writeHandle1(wb, 42);

    WriteDataHandle<int> writeHandle2(&dummyElement, "test2");
    writeHandle2.initialize("same_key");
    BOOST_CHECK_THROW(writeHandle2(wb, 43), std::invalid_argument);
  }

  BOOST_TEST_CHECKPOINT("Test write handle with different type");
  {
    WriteDataHandle<int> writeHandle1(&dummyElement, "test1");
    writeHandle1.initialize("diff_key");
    writeHandle1(wb, 42);

    WriteDataHandle<std::string> writeHandle2(&dummyElement, "test2");
    writeHandle2.initialize("diff_key");
    BOOST_CHECK_THROW(writeHandle2(wb, "test"), std::invalid_argument);
  }

  BOOST_TEST_CHECKPOINT("Test read handle with same type");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("read_key");
    writeHandle(wb, 42);

    ReadDataHandle<int> readHandle1(&dummyElement, "test1");
    readHandle1.initialize("read_key");
    BOOST_CHECK_EQUAL(readHandle1(wb), 42);

    ReadDataHandle<int> readHandle2(&dummyElement, "test2");
    readHandle2.initialize("read_key");
    BOOST_CHECK_EQUAL(readHandle2(wb), 42);
  }

  BOOST_TEST_CHECKPOINT("Test read handle with different type");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("type_key");
    writeHandle(wb, 42);

    ReadDataHandle<std::string> readHandle(&dummyElement, "test");
    readHandle.initialize("type_key");
    BOOST_CHECK_THROW(readHandle(wb), std::out_of_range);
  }
}

BOOST_AUTO_TEST_CASE(WhiteBoardCopy) {
  WhiteBoard wb1;
  WhiteBoard wb2;
  DummySequenceElement dummyElement;

  BOOST_TEST_CHECKPOINT("Test copying from another whiteboard");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("copy_key");
    writeHandle(wb1, 42);

    wb2.copyFrom(wb1);

    ReadDataHandle<int> readHandle(&dummyElement, "test");
    readHandle.initialize("copy_key");
    BOOST_CHECK_EQUAL(readHandle(wb2), 42);
  }

  BOOST_TEST_CHECKPOINT("Test copying with duplicate keys");
  {
    WriteDataHandle<int> writeHandle1(&dummyElement, "test1");
    writeHandle1.initialize("duplicate_key");
    writeHandle1(wb1, 42);

    WriteDataHandle<int> writeHandle2(&dummyElement, "test2");
    writeHandle2.initialize("duplicate_key");
    writeHandle2(wb2, 43);

    BOOST_CHECK_THROW(wb2.copyFrom(wb1), std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(EmulateStateConsistency) {
  DummySequenceElement dummyElement;
  DataHandleBase::StateMapType state;
  WhiteBoard::AliasMapType aliases;
  std::vector<std::unique_ptr<WriteDataHandleBase>> writeHandles;
  std::vector<std::unique_ptr<ReadDataHandleBase>> readHandles;

  BOOST_TEST_CHECKPOINT("Test write handle emulation");
  {
    auto& writeHandle = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<int>>(&dummyElement, "test"));
    writeHandle.initialize("test_key");
    writeHandle.emulate(state, aliases, logger());

    // Verify state after write handle emulation
    BOOST_CHECK(state.contains("test_key"));
    BOOST_CHECK_EQUAL(state["test_key"], &writeHandle);
  }

  BOOST_TEST_CHECKPOINT(
      "Test read handle emulation with compatible write handle");
  {
    auto& readHandle = *readHandles.emplace_back(
        std::make_unique<ReadDataHandle<int>>(&dummyElement, "test"));
    readHandle.initialize("test_key");
    readHandle.emulate(state, aliases, logger());
    // Should not throw as write handle exists with same type
  }

  BOOST_TEST_CHECKPOINT(
      "Test read handle emulation with incompatible write handle");
  {
    state.clear();
    aliases.clear();

    auto& writeHandle = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<std::string>>(&dummyElement, "test"));
    writeHandle.initialize("test_key");
    writeHandle.emulate(state, aliases, logger());

    auto& readHandle = *readHandles.emplace_back(
        std::make_unique<ReadDataHandle<int>>(&dummyElement, "test"));
    readHandle.initialize("test_key");
    ScopedFailureThreshold st(Acts::Logging::Level::FATAL);
    BOOST_CHECK_THROW(readHandle.emulate(state, aliases, logger()),
                      SequenceConfigurationException);
  }

  BOOST_TEST_CHECKPOINT("Test read handle emulation with missing write handle");
  {
    state.clear();
    aliases.clear();
    auto& readHandle = *readHandles.emplace_back(
        std::make_unique<ReadDataHandle<int>>(&dummyElement, "test"));
    readHandle.initialize("missing_key");
    ScopedFailureThreshold st(Acts::Logging::Level::FATAL);
    BOOST_CHECK_THROW(readHandle.emulate(state, aliases, logger()),
                      SequenceConfigurationException);
  }

  BOOST_TEST_CHECKPOINT("Test write handle emulation with duplicate key");
  {
    state.clear();
    aliases.clear();
    auto& writeHandle1 = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<int>>(&dummyElement, "test1"));
    writeHandle1.initialize("duplicate_key");
    writeHandle1.emulate(state, aliases, logger());

    auto& writeHandle2 = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<int>>(&dummyElement, "test2"));
    writeHandle2.initialize("duplicate_key");
    ScopedFailureThreshold st(Acts::Logging::Level::FATAL);
    BOOST_CHECK_THROW(writeHandle2.emulate(state, aliases, logger()),
                      SequenceConfigurationException);
  }

  BOOST_TEST_CHECKPOINT("Test alias handling");
  {
    state.clear();
    aliases.clear();

    auto& writeHandle = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<int>>(&dummyElement, "test"));
    writeHandle.initialize("original_key");
    writeHandle.emulate(state, aliases, logger());

    // Add alias
    aliases.insert({"original_key", "alias_key"});
    state.insert({"alias_key", &writeHandle});

    // Verify read handle works with alias
    auto& readHandle = *readHandles.emplace_back(
        std::make_unique<ReadDataHandle<int>>(&dummyElement, "test"));
    readHandle.initialize("alias_key");
    readHandle.emulate(state, aliases, logger());
    // Should not throw as alias exists and points to compatible handle
  }
}

BOOST_AUTO_TEST_CASE(ConsumeDataHandleTest) {
  WhiteBoard wb;
  DummySequenceElement dummyElement;

  BOOST_TEST_CHECKPOINT("Test basic consume functionality");
  {
    WriteDataHandle<int> writeHandle(&dummyElement, "test");
    writeHandle.initialize("consume_key");
    writeHandle(wb, 42);

    ConsumeDataHandle<int> consumeHandle(&dummyElement, "test");
    consumeHandle.initialize("consume_key");
    BOOST_CHECK_EQUAL(consumeHandle(wb), 42);

    // Verify data is removed after consumption
    BOOST_CHECK(!wb.exists("consume_key"));
    ReadDataHandle<int> readHandle(&dummyElement, "test");
    readHandle.initialize("consume_key");
    BOOST_CHECK_THROW(readHandle(wb), std::out_of_range);
  }

  BOOST_TEST_CHECKPOINT("Test consume handle emulation");
  {
    DataHandleBase::StateMapType state;
    WhiteBoard::AliasMapType aliases;
    std::vector<std::unique_ptr<WriteDataHandleBase>> writeHandles;
    std::vector<std::unique_ptr<ReadDataHandleBase>> readHandles;

    // Add write handle to state
    auto& writeHandle = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<int>>(&dummyElement, "test"));
    writeHandle.initialize("consume_key");
    writeHandle.emulate(state, aliases, logger());

    // Verify consume handle removes key from state
    auto& consumeHandle = *readHandles.emplace_back(
        std::make_unique<ConsumeDataHandle<int>>(&dummyElement, "test"));
    consumeHandle.initialize("consume_key");
    consumeHandle.emulate(state, aliases, logger());
    BOOST_CHECK(!state.contains("consume_key"));

    // Verify another consume handle fails
    auto& consumeHandle2 = *readHandles.emplace_back(
        std::make_unique<ConsumeDataHandle<int>>(&dummyElement, "test2"));
    consumeHandle2.initialize("consume_key");
    ScopedFailureThreshold st(Acts::Logging::Level::FATAL);
    BOOST_CHECK_THROW(consumeHandle2.emulate(state, aliases, logger()),
                      SequenceConfigurationException);
  }

  BOOST_TEST_CHECKPOINT("Test consume handle with incompatible type");
  {
    WriteDataHandle<std::string> writeHandle(&dummyElement, "test");
    writeHandle.initialize("type_key");
    writeHandle(wb, "test string");

    ConsumeDataHandle<int> consumeHandle(&dummyElement, "test");
    consumeHandle.initialize("type_key");
    BOOST_CHECK_THROW(consumeHandle(wb), std::out_of_range);
  }

  BOOST_TEST_CHECKPOINT("Test consume handle with missing data");
  {
    ConsumeDataHandle<int> consumeHandle(&dummyElement, "test");
    consumeHandle.initialize("missing_key");
    BOOST_CHECK_THROW(consumeHandle(wb), std::out_of_range);
  }

  BOOST_TEST_CHECKPOINT("Test consume handle with uninitialized key");
  {
    ConsumeDataHandle<int> consumeHandle(&dummyElement, "test");
    BOOST_CHECK_THROW(consumeHandle(wb), std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(ConsumeDataHandleWithAliases) {
  DummySequenceElement dummyElement;
  DataHandleBase::StateMapType state;
  WhiteBoard::AliasMapType aliases;
  std::vector<std::unique_ptr<WriteDataHandleBase>> writeHandles;
  std::vector<std::unique_ptr<ReadDataHandleBase>> readHandles;

  BOOST_TEST_CHECKPOINT("Test consume handle emulation with aliases");
  {
    aliases.insert({"original_key", "alias_key"});
    // Set up write handle with original key
    auto& writeHandle = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<int>>(&dummyElement, "test"));
    writeHandle.initialize("original_key");
    writeHandle.emulate(state, aliases, logger());

    // Verify initial state
    BOOST_CHECK(state.contains("original_key"));
    BOOST_CHECK(state.contains("alias_key"));
    BOOST_CHECK_EQUAL(state["original_key"], &writeHandle);
    BOOST_CHECK_EQUAL(state["alias_key"], &writeHandle);

    // Emulate consume handle with alias
    auto& consumeHandle = *readHandles.emplace_back(
        std::make_unique<ConsumeDataHandle<int>>(&dummyElement, "test"));
    consumeHandle.initialize("alias_key");
    consumeHandle.emulate(state, aliases, logger());

    // Verify both original and alias keys are removed from state
    BOOST_CHECK(!state.contains("original_key"));
    BOOST_CHECK(!state.contains("alias_key"));
  }

  BOOST_TEST_CHECKPOINT("Test consume handle emulation with original key");
  {
    state.clear();
    aliases.clear();
    aliases.insert({"original_key", "alias_key"});

    // Set up write handle with original key
    auto& writeHandle = *writeHandles.emplace_back(
        std::make_unique<WriteDataHandle<int>>(&dummyElement, "test"));
    writeHandle.initialize("original_key");
    writeHandle.emulate(state, aliases, logger());

    // Verify initial state
    BOOST_CHECK(state.contains("original_key"));
    BOOST_CHECK(state.contains("alias_key"));
    BOOST_CHECK_EQUAL(state["original_key"], &writeHandle);
    BOOST_CHECK_EQUAL(state["alias_key"], &writeHandle);

    // Emulate consume handle with original key
    auto& consumeHandle = *readHandles.emplace_back(
        std::make_unique<ConsumeDataHandle<int>>(&dummyElement, "test"));
    consumeHandle.initialize("original_key");
    consumeHandle.emulate(state, aliases, logger());

    // Verify both original and alias keys are removed from state
    BOOST_CHECK(!state.contains("original_key"));
    BOOST_CHECK(!state.contains("alias_key"));
  }
}

// Custom type with destructor counter for testing
struct DestructorCounter {
  static int count;
  int value;
  explicit DestructorCounter(int v) : value(v) {}
  ~DestructorCounter() { count++; }
};
int DestructorCounter::count = 0;

BOOST_AUTO_TEST_CASE(ConsumeDataHandleDestructor) {
  WhiteBoard wb;
  DummySequenceElement dummyElement;

  BOOST_TEST_CHECKPOINT("Test value destructor is not called when popping");
  {
    // Reset counter
    DestructorCounter::count = 0;

    // Write value to store
    WriteDataHandle<std::unique_ptr<DestructorCounter>> writeHandle(
        &dummyElement, "test");
    writeHandle.initialize("destructor_key");
    writeHandle(wb, std::make_unique<DestructorCounter>(42));

    // Verify initial state
    BOOST_CHECK_EQUAL(DestructorCounter::count, 0);

    // Consume value
    ConsumeDataHandle<std::unique_ptr<DestructorCounter>> consumeHandle(
        &dummyElement, "test");
    consumeHandle.initialize("destructor_key");
    auto value = consumeHandle(wb);

    // Verify value was moved correctly
    BOOST_CHECK_EQUAL(value->value, 42);
    // Verify destructor was not called during pop
    BOOST_CHECK_EQUAL(DestructorCounter::count, 0);

    // Verify value is removed from store
    BOOST_CHECK(!wb.exists("destructor_key"));

    // Value destructor will be called when unique_ptr is destroyed
    value.reset();
    BOOST_CHECK_EQUAL(DestructorCounter::count, 1);
  }
}

BOOST_AUTO_TEST_SUITE_END()
