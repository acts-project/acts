// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/DataHandle.hpp"

namespace ActsExamples {

void WriteDataHandleBase::initialize(const std::string& key) {
  if (key.empty()) {
    throw std::invalid_argument{"Write handle '" + fullName() +
                                "' cannot receive empty key"};
  }
  m_key = key;
}

bool WriteDataHandleBase::isCompatible(const DataHandleBase& other) const {
  return dynamic_cast<const ReadDataHandleBase*>(&other) != nullptr &&
         typeInfo() == other.typeInfo();
}

void ReadDataHandleBase::initialize(const std::string& key) {
  if (key.empty()) {
    throw std::invalid_argument{"Read handle '" + fullName() +
                                "' cannot receive empty key"};
  }
  m_key = key;
}

bool ReadDataHandleBase::isCompatible(const DataHandleBase& other) const {
  return dynamic_cast<const WriteDataHandleBase*>(&other) != nullptr &&
         typeInfo() == other.typeInfo();
}

}  // namespace ActsExamples
