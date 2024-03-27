// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace Acts {

class Volume;
class TrackingVolume;

class BlueprintNode {
 public:
  BlueprintNode(const std::string& name) : m_name(name) {}

  virtual ~BlueprintNode() = default;

  std::string name() const { return m_name; }
  void setName(const std::string& name) { m_name = name; }

  virtual void toStream(std::ostream& os) const;

  virtual Volume& build() = 0;

  // @TODO: This should return the portal "shell"
  virtual void connect(TrackingVolume& parent) = 0;

  void addChild(std::unique_ptr<BlueprintNode> child) {
    m_children.push_back(std::move(child));
  }

  const std::vector<std::unique_ptr<BlueprintNode>>& children() {
    return m_children;
  }

 private:
  std::string m_name;

  std::vector<std::unique_ptr<BlueprintNode>> m_children{};
};

};  // namespace Acts
