// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Mille/MillePedeError.hpp"
#include "ActsPlugins/Mille/MillePedeResultReader.hpp"

#include <format>
#include <fstream>

using namespace ActsPlugins::ActsToMille;

BOOST_AUTO_TEST_SUITE(MillePedeResultReaderSimpleFile)

bool operator==(const mpParameterResult &r1, const mpParameterResult &r2) {
  return (r1.label == r2.label && r1.val == r2.val && r1.start == r2.start &&
          r1.delta == r2.delta && r1.sigma == r2.sigma &&
          r1.nRecords == r2.nRecords);
}

std::ostream &operator<<(std::ostream &str, const mpParameterResult &p) {
  str << std::format(
      " [ label {}, val {}, start {}, delta {}, sigma {}, nRecords {} ]",
      p.label, p.val, p.start, p.delta, p.sigma, p.nRecords);
  return str;
}

/// catch a missing steering file
BOOST_AUTO_TEST_CASE(ReadSimpleFile) {
  std::filesystem::path dummyFileName("testResult.txt");
  std::ofstream dummyResultFile(dummyFileName);
  dummyResultFile << " # this is a comment line and will be ignored."
                  << std::endl;

  std::vector<mpParameterResult> dummyResults{
      {1, 0.01, 0.00, 0.01, 0.001, 144},
      {41, -0.07, 0.00, -0.07, 0.001, 86},
  };
  for (const auto &[label, val, start, delta, sigma, nRec] : dummyResults) {
    dummyResultFile << label << " " << val << " " << start << " " << delta
                    << " " << sigma << " " << nRec << std::endl;
  }
  dummyResultFile.close();
  MillePedeResultReader reader;
  auto res = reader.readParameters(dummyFileName);
  BOOST_CHECK(res.ok());
  BOOST_CHECK_EQUAL(res->size(), dummyResults.size());
  for (std::size_t i = 0; i < res->size(); ++i) {
    BOOST_CHECK(res->at(i) == dummyResults[i]);
  }
}
/// catch a missing steering file
BOOST_AUTO_TEST_CASE(ReadInvalidFile) {
  MillePedeResultReader reader;
  auto res = reader.readParameters("/this/does/hopefully/not/exist?");
  BOOST_CHECK(!res.ok());
  BOOST_CHECK_EQUAL(res.error(), MillePedeError::SolutionNotReadable);
}

BOOST_AUTO_TEST_SUITE_END()
