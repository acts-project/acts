// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/CollectionBaseWriteHandle.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <memory>

#include <podio/CollectionBase.h>

namespace ActsExamples {

CollectionBaseWriteHandle::CollectionBaseWriteHandle(SequenceElement* parent,
                                                     const std::string& name)
    : WriteDataHandleBase{parent, name} {
  registerAsWriteHandle();
}

void CollectionBaseWriteHandle::store(
    WhiteBoard& wb, std::unique_ptr<podio::CollectionBase> collection) const {
  if (!isInitialized()) {
    throw std::runtime_error{"CollectionBaseHandle '" + fullName() +
                             "' not initialized"};
  }

  add(wb, std::move(collection));
}

const std::type_info& CollectionBaseWriteHandle::typeInfo() const {
  return typeid(std::unique_ptr<podio::CollectionBase>);
}

}  // namespace ActsExamples
