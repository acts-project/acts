// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/LinearizedTrack.hpp"

Acts::LinearizedTrack::LinearizedTrack()
{
  m_paramsAtPCA.setZero();
  m_parCovarianceAtPCA.setZero();
  m_linPoint.setZero();
  m_positionJacobian.setZero();
  m_momentumJacobian.setZero();
  m_positionAtPCA.setZero();
  m_momentumAtPCA.setZero();
  m_constTerm.setZero();
}

Acts::LinearizedTrack::LinearizedTrack(
    const Acts::ActsVectorD<5>&    paramsAtPCA,
    const Acts::ActsSymMatrixD<5>& parCovarianceAtPCA,
    const Acts::Vector3D&          linPoint,
    const Acts::ActsMatrixD<5, 3>& positionJacobian,
    const Acts::ActsMatrixD<5, 3>& momentumJacobian,
    const Acts::Vector3D&       positionAtPCA,
    const Acts::Vector3D&       momentumAtPCA,
    const Acts::ActsVectorD<5>& constTerm)
  : m_paramsAtPCA(paramsAtPCA)
  , m_parCovarianceAtPCA(parCovarianceAtPCA)
  , m_linPoint(linPoint)
  , m_positionJacobian(positionJacobian)
  , m_momentumJacobian(momentumJacobian)
  , m_positionAtPCA(positionAtPCA)
  , m_momentumAtPCA(momentumAtPCA)
  , m_constTerm(constTerm)
{
}

Acts::LinearizedTrack&
Acts::LinearizedTrack::operator=(const Acts::LinearizedTrack& other)
{
  if (this != &other) {
    m_paramsAtPCA        = other.m_paramsAtPCA;
    m_parCovarianceAtPCA = other.m_parCovarianceAtPCA;
    m_linPoint           = other.m_linPoint;
    m_positionJacobian   = other.m_positionJacobian;
    m_momentumJacobian   = other.m_momentumJacobian;
    m_positionAtPCA      = other.m_positionAtPCA;
    m_momentumAtPCA      = other.m_momentumAtPCA;
    m_constTerm          = other.m_constTerm;
  }
  return *this;
}