// small change by Lixin to fit pytorch's API

/*
	Matt Dean - 1422434 - mxd434
	
	Goals implemented:
		- Block scan for arbitrary length small vectors - 'blockscan' function
		- Full scan for arbitrary length large vectors	- 'scan' function
			This function decides whether to perform a small (one block) scan or a full (n-level) scan depending on the length of the input vector
		- BCAO for both scans

	Hardware:
		CPU - Intel Core i5-4670k @ 3.4GHz
		GPU - NVIDIA GeForce GTX 760

	Timings:
		10,000,000 Elements
		  host     : 20749 ms
		  gpu      : 7.860768 ms
		  gpu bcao : 4.304064 ms
		
		For more results please see the comment at the bottom of this file

	Extra work:
		Due to the recursive nature of the full scan it can handle n > 3 levels 
	
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"

#include "prefix_sum.h"

// scan.cuh
void sequential_scan(int* output, int* input, int length);
void blockscan(int *output, int *input, int length, bool bcao);
void scan(int *output, int *input, int length, bool bcao);

void scanLargeDeviceArray(int *output, int *input, int length, bool bcao);
void scanSmallDeviceArray(int *d_out, int *d_in, int length, bool bcao);
void scanLargeEvenDeviceArray(int *output, int *input, int length, bool bcao);


// kernels.cuh
__global__ void prescan_arbitrary(int *output, int *input, int n, int powerOfTwo);
__global__ void prescan_arbitrary_unoptimized(int *output, int *input, int n, int powerOfTwo);

__global__ void prescan_large(int *output, int *input, int n, int* sums);
__global__ void prescan_large_unoptimized(int *output, int *input, int n, int *sums);

__global__ void add(int *output, int length, int *n1);
__global__ void add(int *output, int length, int *n1, int *n2);


// utils.h
void _checkCudaError(const char *message, cudaError_t err, const char *caller);
void printResult(const char* prefix, int result, long nanoseconds);
void printResult(const char* prefix, int result, float milliseconds);

bool isPowerOfTwo(int x);
int nextPowerOfTwo(int x);

long get_nanos();

/*///////////////////////////////////*/
/*            New API                */
/*///////////////////////////////////*/

void PrefixSumCUDA(
    const at::Tensor grid_cnt,
    int num_grids,
    at::Tensor grid_off) {
  
  scan(
    grid_off.contiguous().data_ptr<int>(),
    grid_cnt.contiguous().data_ptr<int>(),
    num_grids,
    true
  );

  return;
}

void PrefixSumCPU(
    const at::Tensor grid_cnt,
    int num_grids,
    at::Tensor grid_off) {

  sequential_scan(
    grid_off.contiguous().data_ptr<int>(),
    grid_cnt.contiguous().data_ptr<int>(),
    num_grids
  );

  return;
}


/*///////////////////////////////////*/
/*            scan.cu                */
/*///////////////////////////////////*/
#define checkCudaError(o, l) _checkCudaError(o, l, __func__)

int THREADS_PER_BLOCK = 512;
int ELEMENTS_PER_BLOCK = THREADS_PER_BLOCK * 2;

void sequential_scan(int* output, int* input, int length) {

	output[0] = 0; // since this is a prescan, not a scan
	for (int j = 1; j < length; ++j)
	{
		output[j] = input[j - 1] + output[j - 1];
	}

	return;
}

void blockscan(int *d_out, int *d_in, int length, bool bcao) {
	int powerOfTwo = nextPowerOfTwo(length);
	if (bcao) {
		prescan_arbitrary<<<1, (length + 1) / 2, 2 * powerOfTwo * sizeof(int)>>>(d_out, d_in, length, powerOfTwo);
	}
	else {
		prescan_arbitrary_unoptimized<<<1, (length + 1) / 2, 2 * powerOfTwo * sizeof(int)>>>(d_out, d_in, length, powerOfTwo);
  }
  
  return;
}

void scan(int *d_out, int *d_in, int length, bool bcao) {
	if (length > ELEMENTS_PER_BLOCK) {
		scanLargeDeviceArray(d_out, d_in, length, bcao);
	}
	else {
		scanSmallDeviceArray(d_out, d_in, length, bcao);
	}

	return;
}


void scanLargeDeviceArray(int *d_out, int *d_in, int length, bool bcao) {
	int remainder = length % (ELEMENTS_PER_BLOCK);
	if (remainder == 0) {
		scanLargeEvenDeviceArray(d_out, d_in, length, bcao);
	}
	else {
		// perform a large scan on a compatible multiple of elements
		int lengthMultiple = length - remainder;
		scanLargeEvenDeviceArray(d_out, d_in, lengthMultiple, bcao);

		// scan the remaining elements and add the (inclusive) last element of the large scan to this
		int *startOfOutputArray = &(d_out[lengthMultiple]);
		scanSmallDeviceArray(startOfOutputArray, &(d_in[lengthMultiple]), remainder, bcao);

		add<<<1, remainder>>>(startOfOutputArray, remainder, &(d_in[lengthMultiple - 1]), &(d_out[lengthMultiple - 1]));
	}
}

void scanSmallDeviceArray(int *d_out, int *d_in, int length, bool bcao) {
	int powerOfTwo = nextPowerOfTwo(length);

	if (bcao) {
		prescan_arbitrary<<<1, (length + 1) / 2, 2 * powerOfTwo * sizeof(int)>>>(d_out, d_in, length, powerOfTwo);
	}
	else {
		prescan_arbitrary_unoptimized<<<1, (length + 1) / 2, 2 * powerOfTwo * sizeof(int)>>>(d_out, d_in, length, powerOfTwo);
	}
}

void scanLargeEvenDeviceArray(int *d_out, int *d_in, int length, bool bcao) {
	const int blocks = length / ELEMENTS_PER_BLOCK;
	const int sharedMemArraySize = ELEMENTS_PER_BLOCK * sizeof(int);

	int *d_sums, *d_incr;
	cudaMalloc((void **)&d_sums, blocks * sizeof(int));
	cudaMalloc((void **)&d_incr, blocks * sizeof(int));

	if (bcao) {
		prescan_large<<<blocks, THREADS_PER_BLOCK, 2 * sharedMemArraySize>>>(d_out, d_in, ELEMENTS_PER_BLOCK, d_sums);
	}
	else {
		prescan_large_unoptimized<<<blocks, THREADS_PER_BLOCK, 2 * sharedMemArraySize>>>(d_out, d_in, ELEMENTS_PER_BLOCK, d_sums);
	}

	const int sumsArrThreadsNeeded = (blocks + 1) / 2;
	if (sumsArrThreadsNeeded > THREADS_PER_BLOCK) {
		// perform a large scan on the sums arr
		scanLargeDeviceArray(d_incr, d_sums, blocks, bcao);
	}
	else {
		// only need one block to scan sums arr so can use small scan
		scanSmallDeviceArray(d_incr, d_sums, blocks, bcao);
	}

	add<<<blocks, ELEMENTS_PER_BLOCK>>>(d_out, ELEMENTS_PER_BLOCK, d_incr);

	cudaFree(d_sums);
	cudaFree(d_incr);
}



/*///////////////////////////////////*/
/*            kernels.cu             */
/*///////////////////////////////////*/
#define SHARED_MEMORY_BANKS 32
#define LOG_MEM_BANKS 5

// There were two BCAO optimisations in the paper - this one is fastest
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_MEM_BANKS)

__global__ void prescan_arbitrary(int *output, int *input, int n, int powerOfTwo)
{
	extern __shared__ int temp[];// allocated on invocation
	int threadID = threadIdx.x;

	int ai = threadID;
	int bi = threadID + (n / 2);
	int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
	int bankOffsetB = CONFLICT_FREE_OFFSET(bi);


	if (threadID < n) {
		temp[ai + bankOffsetA] = input[ai];
		temp[bi + bankOffsetB] = input[bi];
	}
	else {
		temp[ai + bankOffsetA] = 0;
		temp[bi + bankOffsetB] = 0;
	}


	int offset = 1;
	for (int d = powerOfTwo >> 1; d > 0; d >>= 1) // build sum in place up the tree
	{
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);

			temp[bi] += temp[ai];
		}
		offset *= 2;
	}

	if (threadID == 0) {
		temp[powerOfTwo - 1 + CONFLICT_FREE_OFFSET(powerOfTwo - 1)] = 0; // clear the last element
	}

	for (int d = 1; d < powerOfTwo; d *= 2) // traverse down tree & build scan
	{
		offset >>= 1;
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);

			int t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
	}
	__syncthreads();

	if (threadID < n) {
		output[ai] = temp[ai + bankOffsetA];
		output[bi] = temp[bi + bankOffsetB];
	}
}

__global__ void prescan_arbitrary_unoptimized(int *output, int *input, int n, int powerOfTwo) {
	extern __shared__ int temp[];// allocated on invocation
	int threadID = threadIdx.x;

	if (threadID < n) {
		temp[2 * threadID] = input[2 * threadID]; // load input into shared memory
		temp[2 * threadID + 1] = input[2 * threadID + 1];
	}
	else {
		temp[2 * threadID] = 0;
		temp[2 * threadID + 1] = 0;
	}


	int offset = 1;
	for (int d = powerOfTwo >> 1; d > 0; d >>= 1) // build sum in place up the tree
	{
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			temp[bi] += temp[ai];
		}
		offset *= 2;
	}

	if (threadID == 0) { temp[powerOfTwo - 1] = 0; } // clear the last element

	for (int d = 1; d < powerOfTwo; d *= 2) // traverse down tree & build scan
	{
		offset >>= 1;
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			int t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
	}
	__syncthreads();

	if (threadID < n) {
		output[2 * threadID] = temp[2 * threadID]; // write results to device memory
		output[2 * threadID + 1] = temp[2 * threadID + 1];
	}
}


__global__ void prescan_large(int *output, int *input, int n, int *sums) {
	extern __shared__ int temp[];

	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int blockOffset = blockID * n;

	int ai = threadID;
	int bi = threadID + (n / 2);
	int bankOffsetA = CONFLICT_FREE_OFFSET(ai);
	int bankOffsetB = CONFLICT_FREE_OFFSET(bi);
	temp[ai + bankOffsetA] = input[blockOffset + ai];
	temp[bi + bankOffsetB] = input[blockOffset + bi];

	int offset = 1;
	for (int d = n >> 1; d > 0; d >>= 1) // build sum in place up the tree
	{
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);

			temp[bi] += temp[ai];
		}
		offset *= 2;
	}
	__syncthreads();


	if (threadID == 0) {
		sums[blockID] = temp[n - 1 + CONFLICT_FREE_OFFSET(n - 1)];
		temp[n - 1 + CONFLICT_FREE_OFFSET(n - 1)] = 0;
	}

	for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
	{
		offset >>= 1;
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);

			int t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
	}
	__syncthreads();

	output[blockOffset + ai] = temp[ai + bankOffsetA];
	output[blockOffset + bi] = temp[bi + bankOffsetB];
}

__global__ void prescan_large_unoptimized(int *output, int *input, int n, int *sums) {
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int blockOffset = blockID * n;

	extern __shared__ int temp[];
	temp[2 * threadID] = input[blockOffset + (2 * threadID)];
	temp[2 * threadID + 1] = input[blockOffset + (2 * threadID) + 1];

	int offset = 1;
	for (int d = n >> 1; d > 0; d >>= 1) // build sum in place up the tree
	{
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			temp[bi] += temp[ai];
		}
		offset *= 2;
	}
	__syncthreads();


	if (threadID == 0) {
		sums[blockID] = temp[n - 1];
		temp[n - 1] = 0;
	}

	for (int d = 1; d < n; d *= 2) // traverse down tree & build scan
	{
		offset >>= 1;
		__syncthreads();
		if (threadID < d)
		{
			int ai = offset * (2 * threadID + 1) - 1;
			int bi = offset * (2 * threadID + 2) - 1;
			int t = temp[ai];
			temp[ai] = temp[bi];
			temp[bi] += t;
		}
	}
	__syncthreads();

	output[blockOffset + (2 * threadID)] = temp[2 * threadID];
	output[blockOffset + (2 * threadID) + 1] = temp[2 * threadID + 1];
}


__global__ void add(int *output, int length, int *n) {
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int blockOffset = blockID * length;

	output[blockOffset + threadID] += n[blockID];
}

__global__ void add(int *output, int length, int *n1, int *n2) {
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int blockOffset = blockID * length;

	output[blockOffset + threadID] += n1[blockID] + n2[blockID];
}


/*///////////////////////////////////*/
/*            utils.cpp              */
/*///////////////////////////////////*/
void _checkCudaError(const char *message, cudaError_t err, const char *caller) {
	if (err != cudaSuccess) {
		fprintf(stderr, "Error in: %s\n", caller);
		fprintf(stderr, message);
		fprintf(stderr, ": %s\n", cudaGetErrorString(err));
		exit(0);
	}
}

void printResult(const char* prefix, int result, long nanoseconds) {
	printf("  ");
	printf(prefix);
	printf(" : %i in %ld ms \n", result, nanoseconds / 1000);
}

void printResult(const char* prefix, int result, float milliseconds) {
	printf("  ");
	printf(prefix);
	printf(" : %i in %f ms \n", result, milliseconds);
}


// from https://stackoverflow.com/a/3638454
bool isPowerOfTwo(int x) {
	return x && !(x & (x - 1));
}

// from https://stackoverflow.com/a/12506181
int nextPowerOfTwo(int x) {
	int power = 1;
	while (power < x) {
		power *= 2;
	}
	return power;
}


// from https://stackoverflow.com/a/36095407
// Get the current time in nanoseconds
long get_nanos() {
	struct timespec ts;
	timespec_get(&ts, TIME_UTC);
	return (long)ts.tv_sec * 1000000000L + ts.tv_nsec;
}
