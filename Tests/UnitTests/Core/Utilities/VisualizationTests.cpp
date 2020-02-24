// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "Acts/Utilities/IVisualization.hpp"
#include "Acts/Utilities/ObjHelper.hpp"
#include "Acts/Utilities/PlyHelper.hpp"

using boost::test_tools::output_test_stream;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_AUTO_TEST_CASE(construction_test) {
  // this doesn't really test anything, other than conformance to the
  // IVisualization interface
  PlyHelper ply;
  ObjHelper obj;

  IVisualization* vis;
  vis = &ply;
  std::cout << *vis << std::endl;
  vis = &obj;
  std::cout << *vis << std::endl;
}

BOOST_AUTO_TEST_CASE(ply_output_test) {
  PlyHelper ply;
  output_test_stream output;

  ply.vertex({0, 0, 0});

  std::string exp = R"(ply
format ascii 1.0
element vertex 1
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
element face 0
property list uchar int vertex_index
element edge 0
property int vertex1
property int vertex2
property uchar red
property uchar green
property uchar blue
end_header
0 0 0 120 120 120
)";

  output << ply;
  BOOST_CHECK(output.is_equal(exp));

  ply.clear();
  ply.vertex({0, 1, 0});

  exp = R"(ply
format ascii 1.0
element vertex 1
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
element face 0
property list uchar int vertex_index
element edge 0
property int vertex1
property int vertex2
property uchar red
property uchar green
property uchar blue
end_header
0 1 0 120 120 120
)";

  output << ply;
  BOOST_CHECK(output.is_equal(exp));

  ply.clear();
  ply.line({0, 0, 1}, {1, 0, 0});

  output << ply;

  exp = R"(ply
format ascii 1.0
element vertex 2
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
element face 0
property list uchar int vertex_index
element edge 1
property int vertex1
property int vertex2
property uchar red
property uchar green
property uchar blue
end_header
0 0 1 120 120 120
1 0 0 120 120 120
0 1 120 120 120
)";

  BOOST_CHECK(output.is_equal(exp));

  ply.clear();
  ply.face({{1, 0, 0}, {1, 1, 0}, {0, 1, 0}});

  output << ply;

  exp = R"(ply
format ascii 1.0
element vertex 3
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
element face 1
property list uchar int vertex_index
element edge 0
property int vertex1
property int vertex2
property uchar red
property uchar green
property uchar blue
end_header
1 0 0 120 120 120
1 1 0 120 120 120
0 1 0 120 120 120
3 0 1 2
)";

  BOOST_CHECK(output.is_equal(exp));
}

BOOST_AUTO_TEST_CASE(obj_output_test) {
  ObjHelper obj;

  output_test_stream output;

  obj.vertex({1, 0, 0});

  output << obj;

  std::string exp = R"(v 1 0 0
)";

  BOOST_CHECK(output.is_equal(exp));

  obj.clear();
  obj.face({{1, 0, 0}, {1, 1, 0}, {0, 1, 0}});
  output << obj;

  exp = R"(v 1 0 0
v 1 1 0
v 0 1 0
f 1 2 3
)";

  BOOST_CHECK(output.is_equal(exp));
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts
