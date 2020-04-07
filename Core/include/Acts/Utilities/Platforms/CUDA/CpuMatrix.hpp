#ifndef CPUMATRIX
#define CPUMATRIX

#include "Acts/Utilities/Platforms/CUDA/CudaMatrix.cu"

// column-major style Matrix Definition

namespace Acts{

template<typename Var_t>
class CudaMatrix;
  
template<typename Var_t>
class CpuMatrix{
  
public:

  CpuMatrix() = default;
  CpuMatrix(size_t nRows, size_t nCols){
    SetSize(nRows,nCols);
    cudaMallocHost(&fHostPtr, fNRows*fNCols*sizeof(Var_t));
  }

  CpuMatrix(size_t nRows, size_t nCols, CudaMatrix<Var_t>* cuMat){
    SetSize(nRows,nCols);
    cudaMallocHost(&fHostPtr, fNRows*fNCols*sizeof(Var_t));
    cudaMemcpy(fHostPtr, cuMat->Get(0,0), fSize*sizeof(Var_t), cudaMemcpyDeviceToHost);   
  }
  
  ~CpuMatrix(){
    cudaFreeHost(fHostPtr);
  }

  void SetSize(size_t row, size_t col){
    fNRows = row;
    fNCols = col;
    fSize  = fNRows*fNCols; 
  }
  
  size_t GetNCols(){ return fNCols; }
  size_t GetNRows(){ return fNRows; }
  size_t GetSize() { return fSize; }
  
  Var_t* Get(size_t row=0, size_t col=0){
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

  void CopyD2H(Var_t* devPtr, size_t len, size_t offset){
    cudaMemcpy(fHostPtr+offset, devPtr, len*sizeof(Var_t), cudaMemcpyDeviceToHost);
  }

  void CopyD2H(Var_t* devPtr, size_t len, size_t offset, cudaStream_t* stream){
    cudaMemcpyAsync(fHostPtr+offset, devPtr, len*sizeof(Var_t), cudaMemcpyDeviceToHost, *stream);
  }
  
private:
  Var_t* fHostPtr;
  size_t fNCols;
  size_t fNRows;
  size_t fSize;
};
  
}

#endif
