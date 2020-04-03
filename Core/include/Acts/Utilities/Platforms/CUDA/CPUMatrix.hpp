#ifndef CPUMATRIX
#define CPUMATRIX

#include "Acts/Utilities/Platforms/CUDA/CPUArray.hpp"
#include "Acts/Utilities/Platforms/CUDA/CUDAMatrix.cu"

// column-major style Matrix Definition

namespace Acts{

template<typename Var_t>
class CUDAMatrix;
  
template<typename Var_t>
class CPUMatrix{
  
public:

  CPUMatrix() = default;
  CPUMatrix(size_t nRows, size_t nCols){
    fNRows = nRows;
    fNCols = nCols;
    cudaMallocHost(&fHostPtr, fNRows*fNCols*sizeof(Var_t));
  }

  CPUMatrix(size_t nRows, size_t nCols, CUDAMatrix<Var_t>* cuMat){
    fNRows = nRows;
    fNCols = nCols;
    fHostPtr = (cuMat->GetCPUArray(fNRows*fNCols,0,0))->Get();
  }
  
  ~CPUMatrix(){
    cudaFreeHost(fHostPtr);
  }

  size_t GetNCols(){ return fNCols; }
  size_t GetNRows(){ return fNRows; }

  Var_t* GetEl(size_t row=0, size_t col=0){
    int offset=row+col*fNRows;
    return fHostPtr+offset;
  }
  
  void SetEl(size_t row, size_t col, Var_t val){
    int offset=row+col*fNRows;
    fHostPtr[offset] = val;
  }
  
  Var_t* GetColumn(size_t col){
    return fHostPtr+col*fNRows;
  }
  Var_t* GetRow(size_t row){    
    Var_t* ret = new Var_t[fNCols];
    for(int i_c=0; i_c<fNCols; i_c++) ret[i_c] = fHostPtr[row+fNRows*i_c];
    return ret;    
  }

  void SetRow(size_t row, Var_t* input){
    for(size_t i_c=0; i_c<fNCols; i_c++){
      fHostPtr[row+fNRows*i_c]=input[i_c];
    }
  }

  void SetColumn(size_t col, Var_t* input){
    fHostPtr[col*fNRows] = input[0];
  }

  void CopyD2H(Var_t* devPtr, size_t len, size_t offset=0){
    cudaMemcpy(fHostPtr, devPtr+offset, len*sizeof(Var_t), cudaMemcpyDeviceToHost);
  }
  
private:
  Var_t* fHostPtr;
  size_t fNCols;
  size_t fNRows;
};
  
}

#endif
