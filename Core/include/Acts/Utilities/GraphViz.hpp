// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/OstreamFormatter.hpp"

#include <format>
#include <ostream>
#include <vector>

namespace Acts::GraphViz {

enum class Shape {
  Box,
  Polygon,
  Ellipse,
  Oval,
  Circle,
  Point,
  Egg,
  Triangle,
  Plaintext,
  Plain,
  Diamond,
  Trapezium,
  Parallelogram,
  House,
  Pentagon,
  Hexagon,
  Septagon,
  Octagon,
  DoubleCircle,
  DoubleOctagon,
  TripleOctagon,
  InvTriangle,
  InvTrapezium,
  InvHouse,
  Mdiamond,
  Msquare,
  Mcircle,
  Rect,
  Rectangle,
  Square,
  Star,
  None,
  Underline,
  Cylinder,
  Note,
  Tab,
  Folder,
  Box3d,
  Component,
  Promoter,
  Cds,
  Terminator,
  Utr,
  PrimerSite,
  RestrictionSite,
  FivePOverhang,
  ThreePOverhang,
  NOverhang,
  Assembly,
  Signature,
  Insulator,
  Ribosite,
  RNAStab,
  ProteaseSite,
  ProteinStab,
  RPromoter,
  RArrow,
  LArrow,
  LPromoter
};

std::ostream& operator<<(std::ostream& os, const Shape& shape);

enum class Style {
  Filled,
  Invisible,
  Diagonals,
  Rounded,
  Dashed,
  Dotted,
  Solid,
  Bold
};

std::ostream& operator<<(std::ostream& os, const Style& style);

/// GraphViz node descriptor.
struct Node {
  /// Node identifier.
  std::string id;
  /// Node label text.
  std::string label = "";
  /// Node shape.
  Shape shape = Shape::Ellipse;
  /// Node style attributes.
  std::vector<Style> style = {Style::Solid};
};

std::ostream& operator<<(std::ostream& os, const Node& node);

/// GraphViz edge descriptor.
struct Edge {
  /// Source node.
  Node from;
  /// Destination node.
  Node to;
  /// Edge style.
  Style style = Style::Solid;
};

std::ostream& operator<<(std::ostream& os, const Edge& node);

}  // namespace Acts::GraphViz

ACTS_OSTREAM_FORMATTER(Acts::GraphViz::Style);
ACTS_OSTREAM_FORMATTER(Acts::GraphViz::Shape);
ACTS_OSTREAM_FORMATTER(Acts::GraphViz::Node);
ACTS_OSTREAM_FORMATTER(Acts::GraphViz::Edge);
