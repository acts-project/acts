// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Mille/ActsToMille.hpp"

#include <memory>

#include "Mille/MilleFactory.h"

using namespace ActsPlugins::ActsToMille;

BOOST_AUTO_TEST_SUITE(ActsMilleBasicCallTests)

BOOST_AUTO_TEST_CASE(OpenMilleRecord) {
  const std::string fname = "MyMilleRecord.root";
  std::unique_ptr<Mille::MilleRecord> theRecord =
      Mille::spawnMilleRecord(fname);
  BOOST_CHECK_NE(theRecord.get(), nullptr);
  theRecord->addData(0.02, 0.005, std::vector<unsigned int>{0u},
                     std::vector<double>{1.0}, std::vector<int>{42},
                     std::vector<double>{-1.0});
  theRecord->writeRecord();
  theRecord->flushOutputFile();
  theRecord.reset(nullptr);

  // read back in
  auto theReader = Mille::spawnMilleReader(fname);
  BOOST_CHECK(theReader->open(fname));
  BOOST_CHECK_NE(theReader.get(), nullptr);
  Mille::MilleBuffer buffer;
  buffer.resize(100, true);
  const auto& [nRead, stat] = theReader->readToBuffer(buffer);
  BOOST_CHECK_EQUAL(stat,
                    4);  // status 4 = successful read with float derivatives
  BOOST_CHECK_EQUAL(
      nRead, 5);  // 5 entries: header + 1 measurement + 1 error + 2 derivatives
  const auto& [nReadSecondCall, statSecondCall] =
      theReader->readToBuffer(buffer);
  BOOST_CHECK_EQUAL(statSecondCall,
                    0);  // status 0 = end of file (should only have one record)
  BOOST_CHECK_EQUAL(nReadSecondCall, 0);  // expect 0 entries read at EOF
}

BOOST_AUTO_TEST_SUITE_END()
