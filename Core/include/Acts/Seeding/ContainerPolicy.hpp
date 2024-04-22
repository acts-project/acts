// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iterator>
#include <vector>

template <typename T>
class ContainerPolicy {
 public:
  virtual void policyInsert(T value) = 0;
};

template <typename T>
class GenericBackInserter {
  ContainerPolicy<T>& policy;

 public:
  using value_type = T;
  using iterator_category = std::output_iterator_tag;
  GenericBackInserter(ContainerPolicy<T>& p) : policy(p) {}

  GenericBackInserter& operator=(const T& value) {
    policy.policyInsert(value);
    return *this;
  }
};

template <typename T>
class VectorPolicy : public ContainerPolicy<T> {
  std::vector<T>& container;

 public:
  VectorPolicy(std::vector<T>& vec) : container(vec) {}

  void policyInsert(T value) override { container.push_back(value); }
};
