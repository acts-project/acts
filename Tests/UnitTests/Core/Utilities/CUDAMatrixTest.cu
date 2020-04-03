
// Developer Note (Beomki Yeo):
// Including PlatformDef.h makes an error...

#include "Acts/Utilities/Platforms/CUDA/CUDAArray.cu"
#include "Acts/Utilities/Platforms/CUDA/CuUtils.cu"
#include "Acts/Utilities/Definitions.hpp"
#include <assert.h>
#include <boost/test/unit_test.hpp>
#include <cuda_profiler_api.h>

namespace Acts{
namespace Test{

template<typename AFloat, int row, int col>
__global__ void MatrixLoadStore_MMA(const Eigen::Matrix<AFloat,row,col>* input,
				Eigen::Matrix<AFloat,row,col>* output){
  int globalId = blockIdx.x * blockDim.x + threadIdx.x;
  output[globalId] = input[globalId];   
}

template<typename AFloat, int row, int col>
__global__ void MatrixLoadStore_AMA(const Eigen::Matrix<AFloat,row,col>* input,
				Eigen::Matrix<AFloat,row,col>* output){
  
  for (int i=0; i<col; i++){
    output[blockIdx.x](threadIdx.x,i) = input[blockIdx.x](threadIdx.x,i);
  }
}
  
  
__global__ void Vector3FLoadStore(const Acts::Vector3F* input, Acts::Vector3F* output){
  int globalId = blockIdx.x * blockDim.x + threadIdx.x;
  output[globalId] = input[globalId];   
}
  
__global__ void Vector3DLoadStore(const Acts::Vector3D* input, Acts::Vector3D* output){
  int globalId = blockIdx.x * blockDim.x + threadIdx.x;
  output[globalId] = input[globalId];   
}
  
BOOST_AUTO_TEST_SUITE(Utilities)
BOOST_AUTO_TEST_CASE( CUDAOBJ_TEST ){
  
  cudaProfilerStart();
  
  using AFloat = float;

  //------------------------------------------------
  // Test Matrix backend
  //------------------------------------------------
  
  const int vecDim = 16;  // Vector Dimension
  const int nVec   = 128; // Number of vectors
    
  dim3 gridSize;
  dim3 blockSize;
  int  bufSize;

  // For misaligned memory access (MMA)
  // For misaligned memory access (MMA)
  // For misaligned memory access (MMA) 
  gridSize  = dim3(1,1,1);
  blockSize = dim3(nVec,1,1);
  bufSize   = gridSize.x * blockSize.x;
  CPUArray<Eigen::Matrix<AFloat, vecDim, 1>>  iMat_MMA_cpu(bufSize);
  for (int i=0; i< bufSize; i++){
    iMat_MMA_cpu[i] = Eigen::Matrix<AFloat,vecDim,1>::Random();
  }

  CUDAArray<Eigen::Matrix<AFloat, vecDim, 1>> iMat_MMA_cuda(bufSize, iMat_MMA_cpu.Get(), bufSize, 0);
  CUDAArray<Eigen::Matrix<AFloat, vecDim, 1>> oMat_MMA_cuda(bufSize);
 
  MatrixLoadStore_MMA<AFloat, vecDim, 1><<< gridSize, blockSize >>>(iMat_MMA_cuda.Get(),oMat_MMA_cuda.Get());
  gpuErrChk ( cudaGetLastError() );
    
  CPUArray<Eigen::Matrix<AFloat, vecDim, 1>> oMat_MMA_cpu = *(oMat_MMA_cuda.GetCPUArray(bufSize));

  // For aligned memory access (AMA)
  // For aligned memory access (AMA)
  // For aligned memory access (AMA)
  gridSize  = dim3(1,1,1);
  blockSize = dim3(vecDim,1,1);
  bufSize   = gridSize.x * blockSize.x;

  CPUArray<Eigen::Matrix<AFloat, vecDim, nVec>>  iMat_AMA_cpu(bufSize);
  for (int i=0; i< bufSize; i++){
    iMat_AMA_cpu[i] = Eigen::Matrix<AFloat,vecDim,nVec>::Random();
  }

  CUDAArray<Eigen::Matrix<AFloat,vecDim,nVec>> iMat_AMA_cuda(bufSize, iMat_AMA_cpu.Get(), bufSize, 0);
  CUDAArray<Eigen::Matrix<AFloat,vecDim,nVec>> oMat_AMA_cuda(bufSize);

  MatrixLoadStore_AMA<AFloat, vecDim, nVec><<< gridSize, blockSize >>>(iMat_AMA_cuda.Get(),oMat_AMA_cuda.Get());

  gpuErrChk ( cudaGetLastError() );
    
  CPUArray<Eigen::Matrix<AFloat,vecDim,nVec>> oMat_AMA_cpu = *(oMat_AMA_cuda.GetCPUArray(bufSize));
  
  cudaProfilerStop();

  for (int i=0; i<nVec; i++){
    BOOST_REQUIRE( iMat_MMA_cpu[i] == oMat_MMA_cpu[i] );
  }
  BOOST_REQUIRE( iMat_AMA_cpu[0] == oMat_AMA_cpu[0] ); 
}
BOOST_AUTO_TEST_SUITE_END()

}
}
