// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>

namespace Acts {

namespace Test {

/// This is a test function that tests the validity of an obj stream
/// It tests for special characters that are not allowed to be contained
///
/// It checks for:
///
/// * Allowed as tags:
///    mtllib, usemtl
///    v, vn, vt, l, f, g, o, s
///    comment lines (starting with #)
/// * v, vn, vt are checked for 3 entries
/// * f is checked for integer > 0 entries only
/// * f is also check for being a triplet if triMesh is true
/// * l is checked for integer > 0 entries only
///
/// @param tString is the test string
/// @param triMesh ist he test if only triangular surfaces exist
///
/// @return a vector of failure messages
inline static std::vector<std::string> testObjString(const std::string& tString,
                                                     bool triMesh = false) {
  std::vector<std::string> errorStrings;
  const std::string w = "[ Invalid obj : ";

  std::vector<std::string> passes = {"#", "usemtl", "mtllib", "o", "g", "s"};

  auto ss = std::stringstream{tString};
  for (std::string line; std::getline(ss, line, '\n');) {
    auto fnbc = line.find_first_not_of(" ", 0);
    if (fnbc != std::string::npos) {
      auto snbc = line.find_first_of(" ", fnbc);
      std::string stag = line.substr(fnbc, snbc - fnbc);

      // Ignore comment line, usemtl, mtllib
      bool pass(false);
      for (const auto& ptag : passes) {
        if (ptag == stag) {
          pass = true;
          break;
        }
      }
      if (pass) {
        continue;
      }
      // Check if vectors are three-dimensional
      auto tnbc = line.find_first_not_of(" ", snbc);
      std::string body = line.substr(tnbc, line.size() - tnbc);

      // Check if we have triplets
      if (stag.find("v") != std::string::npos or
          (stag == std::string("f") and triMesh)) {
        std::vector<std::string> bodySplit;
        boost::split(bodySplit, body, boost::is_any_of(" "));
        if (bodySplit.size() != 3 and stag != std::string("f")) {
          errorStrings.push_back(w + line + " ] " + stag +
                                 " must only have three attributes!");
        } else if (bodySplit.size() != 3) {
          errorStrings.push_back("[ not a triangular mesh : " + line + " ]");
        }
        continue;
      }
      // Check if face and line only have positive integer numbers > 1
      // or deliminator " ", " /"
      if (stag == std::string("f") or stag == std::string("l")) {
        bool onlyDigits =
            (body.find_first_not_of("0123456789/ ") == std::string::npos);
        if (!onlyDigits) {
          errorStrings.push_back(w + line + " ] " + stag +
                                 " can only have positive integers!");
        }
        std::vector<std::string> bodySplit;
        boost::split(bodySplit, body, boost::is_any_of(" "));
        for (auto& bs : bodySplit) {
          if (bs == "0") {
            errorStrings.push_back(w + line +
                                   " ] vertex with index 0 detected!");
          }
        }
        continue;
      }
      errorStrings.push_back(w + line + " ] invalid syntax!");
    }
  }
  return errorStrings;
}

/// This is a test function that tests the validity of an obj stream
/// It tests for special characters that are not allowed to be contained
///
/// It checks for:
///
///
/// @param tString is the test string
/// @param triMesh ist he test if only triangular surfaces exist
///
/// @return a vector of failure messages
inline static std::vector<std::string> testPlyString(const std::string& tString,
                                                     bool triMesh = false) {
  std::vector<std::string> errorStrings;
  const std::string w = "[ Invalid ply : ";

  std::vector<std::string> hPasses = {"format", "element", "property",
                                      "comment"};

  auto ss = std::stringstream{tString};
  size_t lnumber = 0;
  for (std::string line; std::getline(ss, line, '\n');) {
    if (lnumber == 0 and line != "ply") {
      errorStrings.push_back(w + line + " ] first line has to be 'ply");
    }
  }

  (void)triMesh;

  return errorStrings;
}

}  // namespace Test

}  // namespace Acts
