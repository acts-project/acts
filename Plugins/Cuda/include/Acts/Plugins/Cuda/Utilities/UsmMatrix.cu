// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>
#include "cuda.h"
#include "cuda_runtime.h"
#include "Acts/Plugins/Cuda/Utilities/CpuMatrix.hpp"
#include "CudaUtils.cu"

namespace Acts{

template<typename var_t>
class UsmMatrix{

public:
  UsmMatrix()=delete;
  UsmMatrix(size_t nRows, size_t nCols){
    m_setSize(nRows,nCols);
    ACTS_CUDA_ERROR_CHECK( cudaMallocManaged((var_t**)&m_devPtr, m_nRows*m_nCols*sizeof(var_t)) );
    cudaDeviceSynchronize();
  }
  
  ~UsmMatrix(){
    cudaDeviceSynchronize();
    cudaFree(m_devPtr);
  }

  var_t* get(size_t row=0, size_t col=0){
    int offset = row+col*m_nRows;
    return m_devPtr+offset;
  }

  void set(size_t row, size_t col, var_t val) {
    size_t offset = row + col * m_nRows;
    m_devPtr[offset] = val;
  }
  
  void zeros(){
    memset(m_devPtr, 0, m_size*sizeof(var_t));
  }
  
private:
  var_t* m_devPtr; 
  size_t m_nCols;
  size_t m_nRows;
  size_t m_size;
  
  void m_setSize(size_t row, size_t col){
    m_nRows = row;
    m_nCols = col;
    m_size  = m_nRows*m_nCols; 
  }
};

}

