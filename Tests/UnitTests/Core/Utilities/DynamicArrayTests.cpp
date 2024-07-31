// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Utilities/DynamicArray.hpp>

#include <iostream>
#include <sstream>
using namespace Acts;

BOOST_AUTO_TEST_SUITE(DynamicArrays)
BOOST_AUTO_TEST_CASE(DynamicArray_AssignmentCheck) {
  DynamicArray<int, 4> myArray{2, 3, 4, 5};
  int classicalArray[2][3][4][5];

  /// Test first of all the construction...
  // bool fail = false;
  const std::array<size_t, 4>& stripes = myArray.getStripes();
  BOOST_CHECK_EQUAL(stripes[0], 2);
  BOOST_CHECK_EQUAL(stripes[1], 3);
  BOOST_CHECK_EQUAL(stripes[2], 4);
  BOOST_CHECK_EQUAL(stripes[3], 5);

  /// Dummy assignment
  myArray.assign(2);
  for (int& val : myArray) {
    BOOST_CHECK_EQUAL(val, 2);
  }
  /// Assign each element to a different value
  int value = 1;
  for (unsigned int i = 0; i < stripes[0]; ++i) {
    for (unsigned int j = 0; j < stripes[1]; ++j) {
      for (unsigned int k = 0; k < stripes[2]; ++k) {
        for (unsigned int l = 0; l < stripes[3]; ++l) {
          myArray[i][j][k][l] = value;
          classicalArray[i][j][k][l] = value;
          ++value;
        }
      }
    }
  }
  std::stringstream modern_str, classical_str;
  int* array_ptr = reinterpret_cast<int*>(classicalArray);
  for (unsigned int i = 0; i < myArray.size(); ++i) {
    modern_str << myArray.get_val(i) << ",";
    classical_str << *(array_ptr + i) << ",";
  }

  BOOST_CHECK_EQUAL(modern_str.str(), classical_str.str());
  /// check that the assignment did not lead to any overwrite
  value = 1;
  for (unsigned int i = 0; i < stripes[0]; ++i) {
    for (unsigned int j = 0; j < stripes[1]; ++j) {
      for (unsigned int k = 0; k < stripes[2]; ++k) {
        for (unsigned int l = 0; l < stripes[3]; ++l) {
          BOOST_CHECK_EQUAL(myArray[i][j][k][l], value);

          BOOST_CHECK_EQUAL(myArray.get(i, j, k, l), myArray[i][j][k][l]);

          ++value;
        }
      }
    }
  }
  /// Check whether the mapping of global -> local index is bijective
  for (size_t globIdx = 0; globIdx < myArray.size(); ++globIdx) {
    std::array<size_t, 4> axisRanges = myArray.axisIndices(globIdx);
    size_t backIdx = myArray.index(axisRanges);
    BOOST_CHECK_EQUAL(backIdx, globIdx);
  }
  /// Overwrite the stripes
  myArray.changeStripes(5, 4, 3, 2);
  BOOST_CHECK_EQUAL(stripes[0], 5);
  BOOST_CHECK_EQUAL(stripes[1], 4);
  BOOST_CHECK_EQUAL(stripes[2], 3);
  BOOST_CHECK_EQUAL(stripes[3], 2);

  /// Test the exception ability
  myArray.changeStripes(0, 0, 0, 0);
  BOOST_CHECK(!myArray);
  BOOST_CHECK_THROW(myArray.get(0, 0, 0, 0), std::runtime_error);
}
BOOST_AUTO_TEST_SUITE_END()
