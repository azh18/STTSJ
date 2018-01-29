
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "gpuKernel.h"
#include <assert.h>
#include <stdio.h>
#include "GlobalHeader.h"
#include "sm_35_atomic_functions.h"
#include <cuda.h>
#include "curand.h"
#include "device_functions.h"


#define CUDA_CALL(x) { const cudaError_t a = (x); if (a!= cudaSuccess) { printf("\nCUDA Error: %s(err_num=%d)\n", cudaGetErrorString(a), a); cudaDeviceReset(); assert(0);}}
// every task of GPU has a GPUHausInfoTable
typedef struct GPUHausInfoTable {
	uint32_t latlonIdxP, latlonIdxQ; // the first offset of latlon of the traj 
	uint32_t textMatchIdx; //the first offset to set the middle match result 
	uint32_t textualResultIdx; //the first offset to store textual intersect result for each block
	uint32_t pointNumP, pointNumQ; // # point in each traj
	uint32_t textNumP, textNumQ; // total # word in each traj

	// uint32_t idxInTextIdxP, idxInTextIdxQ; //not need, because is same as latlonIdx// idx of first point in textIdx, from this to derive the range of keywords of each point
}GPUHashInfoTable;

typedef struct Latlon {
	float lat;
	float lon;
}Latlon;

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);
void CUDAwarmUp() {
	CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));
	CUDA_CALL(cudaSetDevice(0));
}

//static __inline__ __device__ float atomicMin(float* address, float val) {
//	unsigned long long int* address_as_ull = (unsigned long long int*)address;
//	unsigned long long int old = *address, assumed;
//	assumed = *address_as_ull;
//	do {
//		assumed = old;
//		old = atomicCAS(address_as_ull, assumed, __float_as_longlong(val > __longlong_as_float(assumed) ?
//			__longlong_as_float(assumed) : val));
//	} while (assumed != old);
//	return __longlong_as_float(old);
//}

/*
GPU function definition.
All functions of GPU are defined here.

*/
__global__ void addKernel(int *c, const int *a, const int *b)
{
	int i = threadIdx.x;
	c[i] = a[i] + b[i];
}

__device__ inline float SDistance(float lat1, float lon1, float lat2, float lon2) {
	float latdelta = lat1 - lat2;
	float londelta = lon1 - lon2;
	return sqrt(latdelta*latdelta + londelta*londelta)/ MAX_DIST;
}

__device__ inline float TDistance(int* word1, int* word2, uint32_t wordNum1, uint32_t wordNum2) {
	int tempWords[MAX_KEYWORD_NUM]; uint32_t intersect_size = 0, union_size = 0;
	for (uint32_t idxW = 0; idxW < wordNum1; idxW++) {
		tempWords[idxW] = word1[idxW];
		union_size++;
	}
	for (uint32_t idxW = 0; idxW < wordNum2; idxW++) {
		bool haveSame = false;
		for (uint32_t idxW1 = 0; idxW1 < wordNum1; idxW1++) {
			if (tempWords[idxW1] == word2[idxW]) {
				// intersect_size++;
				haveSame = true;
				break;
			}
		}
		if (haveSame)
			intersect_size++;
		else
			union_size++;
	}
	if (union_size == 0)
		return 0;
	else // do not divide zero!
		return 1.0 - (float)intersect_size / union_size;
}

__global__ void computeHausdorffDistanceByGPUPointLevel(Latlon* latlonP, int* textP, uint32_t* textIdxPArray,
	Latlon* latlonQ, int* textQ, uint32_t* textIdxQArray, int datasizeP, int datasizeQ, 
	uint32_t *wordNumP, uint32_t *wordNumQ, GPUHausInfoTable* taskInfo, float alpha,
	float *HausdorffResult) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	__shared__ GPUHashInfoTable task;
	__shared__ uint32_t pointIdxP, pointNumP, pointIdxQ, pointNumQ;

	// for each thread, get addr and compute distance
	// if traj length is larger than 16, use for loop to process
	__shared__ float tempDist[THREAD_NUM];
	__shared__ float minimunDist[MAX_TRAJ_LENGTH+2]; // one more is because easy to reduce maximum
	__shared__ float maxDist1;

	if (tID == 0) {
		task = taskInfo[bID];
		pointIdxP = task.latlonIdxP;
		pointNumP = task.pointNumP;
		pointIdxQ = task.latlonIdxQ;
		pointNumQ = task.pointNumQ;
	}
	__syncthreads();
	// GPUHashInfoTable task = taskInfo[bID];
	//uint32_t pointIdxP = task.latlonIdxP;
	//uint32_t pointNumP = task.pointNumP;
	//uint32_t pointIdxQ = task.latlonIdxQ;
	//uint32_t pointNumQ = task.pointNumQ;

	if (tID < pointNumP)
		minimunDist[tID] = 9999.0;
	else if (tID < MAX_TRAJ_LENGTH + 2)
		minimunDist[tID] = -1.0;	
	float latP, latQ, lonP, lonQ;
	uint32_t textStartIdxP, textStartIdxQ;
	uint32_t keywordNumP, keywordNumQ;
	for (int i = 0; i < pointNumP; i +=32) {
		if (tID / 8 + i < pointNumP) {
			latP = latlonP[pointIdxP + tID / 8 + i].lat;
			lonP = latlonP[pointIdxP + tID / 8 + i].lon;
			textStartIdxP = textIdxPArray[pointIdxP + tID / 8 + i];
			keywordNumP = wordNumP[pointIdxP + tID / 8 + i];
		}
		for (int j = 0; j < pointNumQ; j += 8) {
			// tID / 16 + i < pointNumP 条件可否去掉？
			if (tID / 8 + i < pointNumP && tID % 8 + j < pointNumQ) {
				latQ = latlonQ[pointIdxQ + tID % 8 + j].lat;
				lonQ = latlonQ[pointIdxQ + tID % 8 + j].lon;
				textStartIdxQ = textIdxQArray[pointIdxQ + tID % 8 + j];
				keywordNumQ = wordNumQ[pointIdxQ + tID % 8 + j];
				// calculate distance
				float TDist = 0;
				//float TDist = TDistance(&textP[textStartIdxP], &textQ[textStartIdxQ], keywordNumP, keywordNumQ);
				float SDist = SDistance(latQ, lonQ, latP, lonP);
				// dist = S*alpha + T*(1-alpha)
				tempDist[tID] = SDist * alpha + TDist * (1 - alpha);
			}
			// compute minimum, loop (may optimize to reduce)
			// only use one thread, low usage!!!
			__syncthreads();
			float minDist = 99999;
			if (tID / 32 == 0 && (tID % 32) < (pointNumP - i)) {
				for (int k = 0; k + j < pointNumQ && k < 8; k++) {
					if (tempDist[k + tID * 8] < minDist)
						minDist = tempDist[k + tID * 8];
				}
				minimunDist[tID % 32 + i] = (minimunDist[tID % 32 + i] < minDist ? minimunDist[tID % 32 + i] : minDist);
				//if (minimunDist[tID % 16 + i] > minDist)
				//	minimunDist[tID % 16 + i] = minDist;
				// atomicMin(&minimunDist[tID % 8 + i], minDist);
			}
			__syncthreads();
		}
	}
	// 归并找最大值
	for (int i = MAX_TRAJ_LENGTH / 2 + 1; i >= 2; ) {

		if (tID < i) {
			minimunDist[tID] = minimunDist[tID] > minimunDist[tID + i] ? minimunDist[tID] : minimunDist[tID + i];
		}
		if (i <= 2) {
			break;
		}
		else
			i = i = (i >> 1) + 1;
		if (tID == 0)
			minimunDist[0] = minimunDist[0] > minimunDist[1] ? minimunDist[0] : minimunDist[1];
		__syncthreads();
	}

	if (tID == 0)
		maxDist1 = minimunDist[0];


	// 计算q对p的距离
	if (tID < pointNumQ)
		minimunDist[tID] = 9999.0;
	else if (tID < MAX_TRAJ_LENGTH + 2)
		minimunDist[tID] = -1.0;

	for (int i = 0; i < pointNumQ; i += 32) {
		if (tID / 8 + i < pointNumQ) {
			latQ = latlonQ[pointIdxQ + tID / 8 + i].lat;
			lonQ = latlonQ[pointIdxQ + tID / 8 + i].lon;
			textStartIdxQ = textIdxQArray[pointIdxQ + tID / 8 + i];
			keywordNumQ = wordNumQ[pointIdxQ + tID / 8 + i];
		}
		for (int j = 0; j < pointNumP; j += 8) {
			// tID / 16 + i < pointNumQ 条件可否去掉？
			if (tID / 8 + i < pointNumQ && tID % 8 + j < pointNumP) {
				latP = latlonP[pointIdxP + tID % 8 + j].lat;
				lonP = latlonP[pointIdxP + tID % 8 + j].lon;
				textStartIdxP = textIdxPArray[pointIdxP + tID % 8 + j];
				keywordNumP = wordNumP[pointIdxP + tID % 8 + j];
				// calculate distance
				float TDist = 0;
				//float TDist = TDistance(&textQ[textStartIdxQ], &textP[textStartIdxP], keywordNumQ, keywordNumP);
				float SDist = SDistance(latQ, lonQ, latP, lonP);
				// dist = S*alpha + T*(1-alpha)
				tempDist[tID] = SDist * alpha + TDist * (1 - alpha);
			}
			// compute minimum, loop (may optimize to reduce)
			// only use a warp, low usage!!!
			__syncthreads();
			float minDist = 9999999;
			if (tID / 32 == 0 && (tID % 32) < (pointNumQ - i)) {
				for (int k = 0; k + j < pointNumP && k < 8; k++) {
					if (tempDist[k + tID * 8] < minDist)
						minDist = tempDist[k + tID * 8];
				}
				minimunDist[tID % 32 + i] = (minimunDist[tID % 32 + i] < minDist ? minimunDist[tID % 32 + i] : minDist);
				// minimunDist[tID % 8 + i] = minDist;

				// atomicMin(&minimunDist[tID % 8 + i], minDist);
				//if (minimunDist[tID % 16 + i] > minDist)
				//	minimunDist[tID % 16 + i] = minDist;
			}
			__syncthreads();
		}
	}
	// 归并找最大值
	// 大量无用归并
	for (int i = MAX_TRAJ_LENGTH / 2+1; i >= 2; ) {

		if (tID < i) {
			minimunDist[tID] = minimunDist[tID] > minimunDist[tID + i] ? minimunDist[tID] : minimunDist[tID + i];
		}
		if (i <= 2) {
			break;
		}
		else
			i = i = (i >> 1) + 1;
		if (tID == 0)
			minimunDist[0] = minimunDist[0] > minimunDist[1] ? minimunDist[0] : minimunDist[1];
		__syncthreads();
	}
	if (tID == 0)
	{
		HausdorffResult[bID] = (minimunDist[0] > maxDist1 ? minimunDist[0] : maxDist1);
	}
	return;
}

/*-----------------------------------------------------------------------------------
	Date: 17/12/20
	keyword-level parallelism version
	solve the problem that some threads are idle when the number of points is small
	use column-oriented like storage to avoid non-coalesced access
	use shared memory to store reduced intermediate results
	leave some space to guarantee no bank conflict
	---------------------------------------------------------------------------------
	waiting for building...
	---------------------------------------------------------------------------------
	matchResult: store data of stage 1 and stage 2, word to word and word to point
	intersectResult: store data of stage 3: point to point
	some data access is not aligned, but just try the performance
	if too many registers are used, divide this kernel function into two functions
*/
__global__ void computeTextualDistanceGPU(Latlon* latlonP, int* textP, uint32_t* textIdxPArray,
	Latlon* latlonQ, int* textQ, uint32_t* textIdxQArray, int datasizeP, int datasizeQ,
	uint32_t *wordNumP, uint32_t *wordNumQ, GPUHausInfoTable* taskInfo, bool* matchResult, float alpha,
	int* intersectResult) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	GPUHashInfoTable task = taskInfo[bID];
	uint32_t pointIdxP = task.latlonIdxP;
	uint32_t pointNumP = task.pointNumP;
	uint32_t pointIdxQ = task.latlonIdxQ;
	uint32_t pointNumQ = task.pointNumQ;
	uint32_t textNumP = task.textNumP;
	uint32_t textNumQ = task.textNumQ;
	uint32_t matchResultIdx = task.textMatchIdx;
	uint32_t textualResultIdx = task.textualResultIdx;

	uint32_t textStartIdxP, textStartIdxQ;
	// parallel keyword match test
	// stage1: pkm&qkn
	for (int i = 0; i < textNumP; i += 8) {
		textStartIdxP = textIdxPArray[pointIdxP] + tID / 32 + i;
		if (textStartIdxP < textNumP) {
			for (int j = 0; j < textNumQ; j += 32) {
				textStartIdxQ = textIdxQArray[pointIdxQ] + tID % 32 + j;
				if (textStartIdxQ<textNumQ)
					matchResult[matchResultIdx + ((textNumQ / 32 + 1) * 32)*(i + tID / 32) + (tID % 32 + j)] = (textP[textStartIdxP] == textQ[textStartIdxQ]);
			}
		}
	}
	__syncthreads();
	// stage2: pkm&q
	int batchSize = 0;
	if (textNumQ <= 32)
		batchSize = 32;
	else if (textNumQ <= 64)
		batchSize = 64;
	else if (textNumQ <= 128)
		batchSize = 128;
	else
		batchSize = 256;
	for (int i = 0; i < textNumQ; i += batchSize) {
		for (int j = 0; j < pointNumP; j += (256 / batchSize)) {
			// add freq in loop
			int pIdx = j + tID / (256 / batchSize); // the idx of point in this traj
			if (pIdx < pointNumP && (i + tID%batchSize < textNumQ)) { // filter threads
				int totalMatch = 0;
				for (int k = 0; k < wordNumP[pIdx + pointIdxP]; k++) {
					totalMatch += matchResult[matchResultIdx + (textNumQ / 32 + 1) * 32 * (textIdxPArray[pIdx + pointIdxP] - textIdxPArray[pointIdxP] + k) + i + tID%batchSize];
				}
				// store final match num in the first line of matchResult matrix
				matchResult[matchResultIdx + (textNumQ / 32 + 1) * 32 * (textIdxPArray[pIdx + pointIdxP] - textIdxPArray[pointIdxP]) + i + tID%batchSize] = totalMatch;
			}
		}
	}
	__syncthreads();
	// stage3:p&q
	for (int i = 0; i < pointNumP; i += 32) {
		for (int j = 0; j < pointNumQ; j += 8) {
			if (i + tID % 32 < pointNumP && j + tID / 32 < pointNumQ) {
				int pIdx = i + tID % 32;
				int qIdx = j + tID / 32;
				int iterNum = wordNumQ[qIdx + pointIdxQ];
				int intersectNum = 0;
				for (int k = 0; k < iterNum; k++) {
					intersectNum += matchResult[matchResultIdx + (textNumQ / 32 + 1) * 32 * (textIdxPArray[pIdx + pointIdxP] - textIdxPArray[pointIdxP]) + k + textIdxQArray[qIdx + pointIdxQ] - textIdxQArray[pointIdxQ]];
				}
				intersectResult[textualResultIdx + pointNumQ*pIdx + qIdx] = intersectNum;
			}
		}
	}
	return;
}

__global__ void computeHausdorffDistanceByGPUKeywordLevel(Latlon* latlonP, int* textP, uint32_t* textIdxPArray,
	Latlon* latlonQ, int* textQ, uint32_t* textIdxQArray, int datasizeP, int datasizeQ,
	uint32_t *wordNumP, uint32_t *wordNumQ, GPUHausInfoTable* taskInfo, bool* matchResult, float alpha,
	float *HausdorffResult, int* intersectResult) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	GPUHashInfoTable task = taskInfo[bID];
	uint32_t pointIdxP = task.latlonIdxP;
	uint32_t pointNumP = task.pointNumP;
	uint32_t pointIdxQ = task.latlonIdxQ;
	uint32_t pointNumQ = task.pointNumQ;
	uint32_t textNumP = task.textNumP;
	uint32_t textNumQ = task.textNumQ;
	uint32_t matchResultIdx = task.textMatchIdx;

	uint32_t textStartIdxP, textStartIdxQ;
	// parallel keyword match test
	// stage1: pkm&qkn
	for (int i = 0; i < textNumP; i += 8) {
		textStartIdxP = textIdxPArray[pointIdxP] + tID / 32 + i;
		if (textStartIdxP < textNumP) {
			for (int j = 0; j < textNumQ; j += 32) {
				textStartIdxQ = textIdxQArray[pointIdxQ] + tID % 32 + j;
				if(textStartIdxQ<textNumQ)
					matchResult[matchResultIdx + ((textNumQ / 32 + 1) * 32)*(i + tID / 32) + (tID % 32 + j)] = (textP[textStartIdxP] == textQ[textStartIdxQ]);
			}
		}
	}
	__syncthreads();
	// stage2: pkm&q
	int batchSize = 0;
	if (textNumQ <= 32)
		batchSize = 32;
	else if (textNumQ <= 64)
		batchSize = 64;
	else if (textNumQ <= 128)
		batchSize = 128;
	else
		batchSize = 256;
	for (int i = 0; i < textNumQ; i += batchSize) {
		for (int j = 0; j < pointNumP; j += (256 / batchSize)) {
			// add freq in loop
			int pIdx = j + tID / (256 / batchSize); // the idx of point in this traj
			if (pIdx < pointNumP && (i + tID%batchSize < textNumQ)) { // filter threads
				int totalMatch = 0;
				for (int k = 0; k < wordNumP[pIdx + pointIdxP]; k++) {
					totalMatch += matchResult[matchResultIdx + (textNumQ / 32 + 1) * 32 * (textIdxPArray[pIdx + pointIdxP] - textIdxPArray[pointIdxP] + k) + i + tID%batchSize];
				}
				// store final match num in the first line of matchResult matrix
				matchResult[matchResultIdx + (textNumQ / 32 + 1) * 32 * (textIdxPArray[pIdx + pointIdxP] - textIdxPArray[pointIdxP]) + i + tID%batchSize] = totalMatch;
			}
		}
	}
	__syncthreads();
	// stage3:p&q
	for (int i = 0; i < pointNumP; i += 32) {
		for (int j = 0; j < pointNumQ; j += 8) {
			if (i + tID % 32 < pointNumP && j + tID/32 < pointNumQ) {
				int pIdx = i + tID % 32;
				int qIdx = j + tID / 32;
				int iterNum = wordNumQ[qIdx + pointIdxQ];
				int intersectNum = 0;
				for (int k = 0; k < iterNum; k++) {
					intersectNum += matchResult[matchResultIdx + (textNumQ / 32 + 1) * 32 * (textIdxPArray[pIdx + pointIdxP] - textIdxPArray[pointIdxP]) + k + textIdxQArray[qIdx + pointIdxQ] - textIdxQArray[pointIdxQ]];
				}
				intersectResult[pointNumQ*pIdx + qIdx] = intersectNum;
			}
		}
	}
	// outputresult





	// for each thread, get addr and compute distance
	// if traj length is larger than 16, use for loop to process
	__shared__ float tempDist[THREAD_NUM];
	__shared__ float minimunDist[MAX_TRAJ_LENGTH + 2]; // one more is because easy to reduce maximum
	__shared__ float maxDist1;

	if (tID < pointNumP)
		minimunDist[tID] = 9999.0;
	else if (tID < MAX_TRAJ_LENGTH + 2)
		minimunDist[tID] = -1.0;
	float latP, latQ, lonP, lonQ;
	uint32_t keywordNumP, keywordNumQ;
	for (int i = 0; i < pointNumP; i += 16) {
		if (tID / 16 + i < pointNumP) {
			latP = latlonP[pointIdxP + tID / 16 + i].lat;
			lonP = latlonP[pointIdxP + tID / 16 + i].lon;
			textStartIdxP = textIdxPArray[pointIdxP + tID / 16 + i];
			keywordNumP = wordNumP[pointIdxP + tID / 16 + i];
		}
		for (int j = 0; j < pointNumQ; j += 16) {
			// tID / 16 + i < pointNumP 条件可否去掉？
			if (tID / 16 + i < pointNumP && tID % 16 + j < pointNumQ) {
				latQ = latlonQ[pointIdxQ + tID % 16 + j].lat;
				lonQ = latlonQ[pointIdxQ + tID % 16 + j].lon;
				textStartIdxQ = textIdxQArray[pointIdxQ + tID % 16 + j];
				keywordNumQ = wordNumQ[pointIdxQ + tID % 16 + j];
				// calculate distance
				float TDist = TDistance(&textP[textStartIdxP], &textQ[textStartIdxQ], keywordNumP, keywordNumQ);
				float SDist = SDistance(latQ, lonQ, latP, lonP);
				// dist = S*alpha + T*(1-alpha)
				tempDist[tID] = SDist * alpha + TDist * (1 - alpha);
			}
			// compute minimum, loop (may optimize to reduce)
			// only use a warp
			__syncthreads();
			float minDist = 99999;
			if (tID / 16 == 0 && (tID % 16) < (pointNumP - i)) {
				for (int k = 0; k + j < pointNumQ && k < 16; k++) {
					if (tempDist[k + tID * 16] < minDist)
						minDist = tempDist[k + tID * 16];
				}
				minimunDist[tID % 16 + i] = (minimunDist[tID % 16 + i] < minDist ? minimunDist[tID % 16 + i] : minDist);
				//if (minimunDist[tID % 16 + i] > minDist)
				//	minimunDist[tID % 16 + i] = minDist;
				// atomicMin(&minimunDist[tID % 16 + i], minDist);
			}
			__syncthreads();
		}
	}
	// 归并找最大值
	for (int i = MAX_TRAJ_LENGTH / 2 + 1; i >= 2; ) {

		if (tID < i) {
			minimunDist[tID] = minimunDist[tID] > minimunDist[tID + i] ? minimunDist[tID] : minimunDist[tID + i];
		}
		if (i <= 2) {
			break;
		}
		else
			i = i = (i >> 1) + 1;
		if (tID == 0)
			minimunDist[0] = minimunDist[0] > minimunDist[1] ? minimunDist[0] : minimunDist[1];
		__syncthreads();
	}

	if (tID == 0)
		maxDist1 = minimunDist[0];


	// 计算q对p的距离
	if (tID < pointNumQ)
		minimunDist[tID] = 9999.0;
	else if (tID < MAX_TRAJ_LENGTH + 2)
		minimunDist[tID] = -1.0;

	for (int i = 0; i < pointNumQ; i += 16) {
		if (tID / 16 + i < pointNumQ) {
			latQ = latlonQ[pointIdxQ + tID / 16 + i].lat;
			lonQ = latlonQ[pointIdxQ + tID / 16 + i].lon;
			textStartIdxQ = textIdxQArray[pointIdxQ + tID / 16 + i];
			keywordNumQ = wordNumQ[pointIdxQ + tID / 16 + i];
		}
		for (int j = 0; j < pointNumP; j += 16) {
			// tID / 16 + i < pointNumQ 条件可否去掉？
			if (tID / 16 + i < pointNumQ && tID % 16 + j < pointNumP) {
				latP = latlonP[pointIdxP + tID % 16 + j].lat;
				lonP = latlonP[pointIdxP + tID % 16 + j].lon;
				textStartIdxP = textIdxPArray[pointIdxP + tID % 16 + j];
				keywordNumP = wordNumP[pointIdxP + tID % 16 + j];
				// calculate distance
				float TDist = TDistance(&textQ[textStartIdxQ], &textP[textStartIdxP], keywordNumQ, keywordNumP);
				float SDist = SDistance(latQ, lonQ, latP, lonP);
				// dist = S*alpha + T*(1-alpha)
				tempDist[tID] = SDist * alpha + TDist * (1 - alpha);
			}
			// compute minimum, loop (may optimize to reduce)
			// only use a warp
			__syncthreads();
			float minDist = 9999999;
			if (tID / 16 == 0 && (tID % 16) < (pointNumQ - i)) {
				for (int k = 0; k + j < pointNumP && k < 16; k++) {
					if (tempDist[k + tID * 16] < minDist)
						minDist = tempDist[k + tID * 16];
				}
				minimunDist[tID % 16 + i] = (minimunDist[tID % 16 + i] < minDist ? minimunDist[tID % 16 + i] : minDist);
				// atomicMin(&minimunDist[tID % 16 + i], minDist);
				//if (minimunDist[tID % 16 + i] > minDist)
				//	minimunDist[tID % 16 + i] = minDist;
			}
			__syncthreads();
		}
	}
	// 归并找最大值
	for (int i = MAX_TRAJ_LENGTH / 2 + 1; i >= 2; ) {

		if (tID < i) {
			minimunDist[tID] = minimunDist[tID] > minimunDist[tID + i] ? minimunDist[tID] : minimunDist[tID + i];
		}
		if (i <= 2) {
			break;
		}
		else
			i = i = (i >> 1) + 1;
		if (tID == 0)
			minimunDist[0] = minimunDist[0] > minimunDist[1] ? minimunDist[0] : minimunDist[1];
		__syncthreads();
	}
	if (tID == 0)
	{
		HausdorffResult[bID] = (minimunDist[0] > maxDist1 ? minimunDist[0] : maxDist1);
	}
	return;
}

/*
CPU function definition.
All functions of CPU are defined here.

*/

size_t calculateDatasize_TrajSet(vector<STTraj> &trajSet) {
	size_t datasize = 0;
	for (vector<STTraj>::iterator it = trajSet.begin(); it != trajSet.end(); it++) {
		for (vector<STPoint>::iterator itp = it->points.begin(); itp != it->points.end(); itp++) {
			datasize += (2 * sizeof(float) + (itp->keywords.size()) * sizeof(int));
		}
	}
	return datasize;
}

// return the bytes copied, from pStart
// numKeywords: #keywords of each point
size_t copySTTrajToArray(STTraj &traj, char* pStart, size_t *numKeywords) {
	vector<STPoint>::iterator itp;
	size_t ptCnt = 0;
	char *s = pStart;
	for (itp = traj.points.begin(); itp != traj.points.end(); itp++) {
		numKeywords[ptCnt++] = itp->keywords.size();
		memcpy(pStart, &itp->lat, sizeof(float));
		pStart += sizeof(float);
		memcpy(pStart, &itp->lon, sizeof(float));
		pStart += sizeof(float);
		for (vector<int>::iterator itk = itp->keywords.begin(); itk != itp->keywords.end(); itk++) {
			memcpy(pStart, &(*itk), sizeof(int));
			pStart += sizeof(int);
		}
	}
	return (pStart - s);
}

void* GPUMalloc(size_t byteNum) {
	void *addr;
	CUDA_CALL(cudaMalloc((void**)&addr, byteNum));
	return addr;
}



int calculateDistanceGPU(vector<STTraj> &trajSetP,
	vector<STTraj> &trajSetQ,
	map<trajPair, float> &result,
	void* baseGPUAddr, void* baseGPUAddr8byteAligned, float alpha, float epsilon,
	cudaStream_t &stream, int parallelLevel) {
	// Latlon *latlonDataPCPU, *latlonDataQCPU; // latlon array
	// int *textDataPCPU, *textDataQCPU; 
	vector<int> textDataPCPU, textDataQCPU; // keyword array
	vector<uint32_t> textIdxPCPU, textIdxQCPU; // keyword idx for each point (to locate the where and how many keywords for a point)
	vector<uint32_t> numWordPCPU, numWordQCPU; // keyword num in each point
	size_t dataSizeP = trajSetP.size(), dataSizeQ = trajSetQ.size();
	// 所有点线性排列
	vector<Latlon> latlonDataPCPU, latlonDataQCPU;
	GPUHausInfoTable *hausTaskInfoCPU = (GPUHashInfoTable*)malloc(sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ);

	MyTimer timer;
	timer.start();
	//Latlon *latlon = latlonDataPCPU;
	//size_t keywordNum[1000];
	uint32_t pointCnt = 0,  textCnt = 0;
	// process set P
	for (size_t i = 0; i < trajSetP.size(); i++) {
		//update table for P[i][]
		uint32_t totalWordNum = 0;
		for (size_t j = i*dataSizeQ; j < (i+1)*dataSizeQ; j++) {
			hausTaskInfoCPU[j].pointNumP = (uint32_t)trajSetP[i].points.size();
			hausTaskInfoCPU[j].latlonIdxP = pointCnt;
		}
		//insert data and update idx
		for (size_t j = 0; j < trajSetP[i].points.size(); j++) {
			Latlon p;
			p.lat = trajSetP[i].points[j].lat;
			p.lon = trajSetP[i].points[j].lon;
			latlonDataPCPU.push_back(p);
			numWordPCPU.push_back((uint32_t)trajSetP[i].points[j].keywords.size());
			pointCnt++;
			textIdxPCPU.push_back(textCnt);
			totalWordNum += (uint32_t)trajSetP[i].points[j].keywords.size();
			for (size_t k = 0; k < trajSetP[i].points[j].keywords.size(); k++) {
				textDataPCPU.push_back(trajSetP[i].points[j].keywords[k]);
				textCnt++; // textCnt is responsible for record textIdx
			}
		}
		for (size_t j = i*dataSizeQ; j < (i + 1)*dataSizeQ; j++) {
			hausTaskInfoCPU[j].textNumP = totalWordNum;
		}
		// when insert keywords, guarantee that the first word idx a traj
		// is aligned with 32 bytes. Push data directly in GPU.
		if (textCnt % 8 != 0) {
			uint32_t cntDiff = (textCnt / 8 + 1) * 8 - textCnt;
			for (uint32_t m = 0; m < cntDiff; m++)
				textDataPCPU.push_back(-1);
			textCnt += cntDiff;
		}
		// printf("textIdx:%d\t", textCnt);
	}
	void* latlonDataPGPU, *latlonDataQGPU, *textDataPGPU, *textDataQGPU, *textIdxPGPU, *textIdxQGPU, *numWordPGPU, *numWordQGPU;
	void *pNow = baseGPUAddr, *pNow8ByteAligned = baseGPUAddr8byteAligned;
	// Copy data of P to GPU
	CUDA_CALL(cudaMemcpyAsync(pNow8ByteAligned, &latlonDataPCPU[0], sizeof(Latlon)*latlonDataPCPU.size(), cudaMemcpyHostToDevice, stream));
	latlonDataPGPU = pNow8ByteAligned;
	pNow8ByteAligned = (void*)((char*)pNow8ByteAligned + (sizeof(Latlon)*latlonDataPCPU.size() / 32 + 1) * 32);
	CUDA_CALL(cudaMemcpyAsync(pNow, &textDataPCPU[0], sizeof(int)*textDataPCPU.size(), cudaMemcpyHostToDevice, stream));
	textDataPGPU = pNow;
	pNow = (void*)((char*)pNow + (sizeof(int)*textDataPCPU.size() / 32 + 1) * 32);
	CUDA_CALL(cudaMemcpyAsync(pNow, &textIdxPCPU[0], sizeof(uint32_t)*textIdxPCPU.size(), cudaMemcpyHostToDevice, stream));
	textIdxPGPU = pNow;
	pNow = (void*)((char*)pNow + (sizeof(uint32_t)*textIdxPCPU.size() / 32 + 1) * 32);
	CUDA_CALL(cudaMemcpyAsync(pNow, &numWordPCPU[0], sizeof(uint32_t)*numWordPCPU.size(), cudaMemcpyHostToDevice, stream));
	numWordPGPU = pNow;
	pNow = (void*)((char*)pNow + (sizeof(uint32_t)*numWordPCPU.size() / 32 + 1) * 32);

	// process set Q
	pointCnt = 0; textCnt = 0;
	for (size_t i = 0; i < trajSetQ.size(); i++) {
		//update table for Q[][i]
		uint32_t totalWordNum = 0;
		for (size_t j = i; j < dataSizeP*dataSizeQ; j+=dataSizeQ) {
			hausTaskInfoCPU[j].pointNumQ = (uint32_t)trajSetQ[i].points.size();
			hausTaskInfoCPU[j].latlonIdxQ = pointCnt;
		}
		//insert data and update idx
		for (size_t j = 0; j < trajSetQ[i].points.size(); j++) {
			Latlon p;
			p.lat = trajSetQ[i].points[j].lat;
			p.lon = trajSetQ[i].points[j].lon;
			// printf("%d,%d,",i, j);
			latlonDataQCPU.push_back(p);
			numWordQCPU.push_back((uint32_t)trajSetQ[i].points[j].keywords.size());
			totalWordNum += (uint32_t)trajSetQ[i].points[j].keywords.size();
			pointCnt++;
			textIdxQCPU.push_back(textCnt);
			for (size_t k = 0; k < trajSetQ[i].points[j].keywords.size(); k++) {
				textDataQCPU.push_back(trajSetQ[i].points[j].keywords[k]);
				textCnt++;
			}
		}
		for (size_t j = i; j < dataSizeP*dataSizeQ; j += dataSizeQ) {
			hausTaskInfoCPU[j].textNumQ = totalWordNum;
		}
		// when insert keywords, guarantee that the first word idx a traj
		// is aligned with 32 bytes. Push data directly in GPU.
		if (textCnt % 8 != 0) {
			uint32_t cntDiff = (textCnt / 8 + 1) * 8 - textCnt;
			for (uint32_t m = 0; m < cntDiff; m++)
				textDataQCPU.push_back(-1);
			textCnt += cntDiff;
		}
	}

	uint32_t matchResultIdx = 0; // the idx of the match matrix for one trajectory pair
	uint32_t textualIntersectResultIdx = 0; // the idx of the intersect result for each block
	for (uint32_t i = 0; i < dataSizeP*dataSizeQ; i++) {
		hausTaskInfoCPU[i].textMatchIdx = matchResultIdx;
		hausTaskInfoCPU[i].textualResultIdx = textualIntersectResultIdx;
		uint32_t textNumP = hausTaskInfoCPU[i].textNumP;
		if (textNumP % 8 != 0) {
			uint32_t cntDiff = (textNumP / 8 + 1) * 8 - textNumP;
			textNumP += cntDiff;
		}
		// guarantee textNumP is the multiple of 8
		matchResultIdx += textNumP*hausTaskInfoCPU[i].textNumQ;
		textualIntersectResultIdx += hausTaskInfoCPU[i].pointNumP*hausTaskInfoCPU[i].pointNumQ;
	}
	// after this loop, textualIntersectResultIdx becomes the total number of elements in textualIntersectResult


	// cpu data build finished
	// transfer data to GPU
	// Copy data of P to GPU
	// this order is to guarantee aligned load in GPU

	CUDA_CALL(cudaMemcpyAsync(pNow8ByteAligned, &latlonDataQCPU[0], sizeof(Latlon)*latlonDataQCPU.size(), cudaMemcpyHostToDevice, stream));
	latlonDataQGPU = pNow8ByteAligned;
	pNow8ByteAligned = (void*)((char*)pNow8ByteAligned + (sizeof(Latlon)*latlonDataQCPU.size() / 32 + 1) * 32);

	CUDA_CALL(cudaMemcpyAsync(pNow, &textDataQCPU[0], sizeof(int)*textDataQCPU.size(), cudaMemcpyHostToDevice, stream));
	textDataQGPU = pNow;
	pNow = (void*)((char*)pNow + (sizeof(int)*textDataQCPU.size() / 32 + 1) * 32);

	CUDA_CALL(cudaMemcpyAsync(pNow, &textIdxQCPU[0], sizeof(uint32_t)*textIdxQCPU.size(), cudaMemcpyHostToDevice, stream));
	textIdxQGPU = pNow;
	pNow = (void*)((char*)pNow + (sizeof(uint32_t)*textIdxQCPU.size() / 32 + 1) * 32);
	CUDA_CALL(cudaMemcpyAsync(pNow, &numWordQCPU[0], sizeof(uint32_t)*numWordQCPU.size(), cudaMemcpyHostToDevice, stream));
	numWordQGPU = pNow;
	pNow = (void*)((char*)pNow + (sizeof(uint32_t)*numWordQCPU.size() / 32 + 1) * 32);

	void *taskinfoTable;
	CUDA_CALL(cudaMemcpyAsync(pNow, hausTaskInfoCPU, sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ, cudaMemcpyHostToDevice, stream));
	taskinfoTable = pNow;
	pNow = (void*)((char*)pNow + (sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ / 32 + 1) * 32);
	
	// reserve memory for hit table: size: numWordP*numWordQ
	void *textMatchResult = pNow;

	// data copy finish, invoke kernel
	float *distanceResult,*distanceResultGPU;
	CUDA_CALL(cudaHostAlloc((void**)&distanceResult, dataSizeP*dataSizeQ*sizeof(float), cudaHostAllocMapped));
	CUDA_CALL(cudaHostGetDevicePointer((void**)&distanceResultGPU, distanceResult, 0));
	int *textualIntersectResult, *textualIntersectResultGPU;
	CUDA_CALL(cudaHostAlloc((void**)&textualIntersectResult, textualIntersectResultIdx * sizeof(int), cudaHostAllocMapped));
	CUDA_CALL(cudaHostGetDevicePointer((void**)&textualIntersectResultGPU, textualIntersectResult, 0));
	//cudaDeviceSynchronize();
	//timer.stop();
	//std::cout << timer.elapse() << "ms to copy data" << std::endl;

	cudaDeviceSynchronize();
	timer.stop();
	printf("stage 1: %f ms\n", timer.elapse());

	// compute textual distance
	//timer.start();
	//computeTextualDistanceGPU << <(uint32_t)dataSizeP*(uint32_t)dataSizeQ, THREAD_NUM, 0, stream >> > ((Latlon*)latlonDataPGPU, (int*)textDataPGPU, (uint32_t*)textIdxPGPU,
	//	(Latlon*)latlonDataQGPU, (int*)textDataQGPU, (uint32_t*)textIdxQGPU, (int)dataSizeP, (int)dataSizeQ,
	//	(uint32_t*)numWordPGPU, (uint32_t*)numWordQGPU, (GPUHashInfoTable*)taskinfoTable, (bool*)textMatchResult, alpha,
	//	(int*)textualIntersectResultGPU);
	//cudaDeviceSynchronize();
	//timer.stop();
	//printf("stage 2: %f ms\n", timer.elapse());
	//for (int i = 0; i < dataSizeP*dataSizeQ; i++) {
	//	uint32_t textualDistIdx = hausTaskInfoCPU[i].textualResultIdx;
	//	int Ptraj = i / dataSizeQ;
	//	int Qtraj = i % dataSizeQ;
	//	printf("\nTDist of T(%d) and T(%d):\n", Ptraj, Qtraj);
	//	for (int k = textualDistIdx; k < textualDistIdx + hausTaskInfoCPU[i].pointNumP*hausTaskInfoCPU[i].pointNumQ; k++) {
	//		printf("%d,", textualIntersectResult[k]);
	//	}
	//}

	timer.start();
	computeHausdorffDistanceByGPUPointLevel << <(uint32_t)dataSizeP*(uint32_t)dataSizeQ, THREAD_NUM, 0, stream >> > ((Latlon*)latlonDataPGPU, 
		(int*)textDataPGPU, (uint32_t*)textIdxPGPU, (Latlon*)latlonDataQGPU, 
		(int*)textDataQGPU, (uint32_t*)textIdxQGPU, 
		(int)dataSizeP, (int)dataSizeQ,
		(uint32_t*)numWordPGPU, (uint32_t*)numWordQGPU, 
		(GPUHashInfoTable*)taskinfoTable, alpha,
		(float*)distanceResultGPU);


	//computeHausdorffDistanceByGPUKeywordLevel << <(uint32_t)dataSizeP*(uint32_t)dataSizeQ, THREAD_NUM, 0, stream >> > ((Latlon*)latlonDataPGPU,
	//	(int*)textDataPGPU, (uint32_t*)textIdxPGPU, (Latlon*)latlonDataQGPU,
	//	(int*)textDataQGPU, (uint32_t*)textIdxQGPU,
	//	(int)dataSizeP, (int)dataSizeQ,
	//	(uint32_t*)numWordPGPU, (uint32_t*)numWordQGPU,
	//	(GPUHashInfoTable*)taskinfoTable, (bool*)textMatchResult, alpha,
	//	(float*)distanceResultGPU);


	cudaDeviceSynchronize();
	timer.stop();
	printf("stage 3: %f ms\n", timer.elapse());
	//for (int i = 0; i < dataSizeP*dataSizeQ; i++) {
	//	printf("d(%zd,%zd)=%f\t", i / dataSizeQ, i%dataSizeQ, distanceResult[i]);
	//}
	// write result
	for (int i = 0; i < dataSizeP*dataSizeQ; i++) {
		if(distanceResult[i] <= epsilon)
			result[trajPair(i / dataSizeQ, i%dataSizeQ)] = distanceResult[i];
	}

	// free memory
	CUDA_CALL(cudaFreeHost(distanceResult));
	free(hausTaskInfoCPU);

	return 0;
}



 
/*
int main()
{
    const int arraySize = 5;
    const int a[arraySize] = { 1, 2, 3, 4, 5 };
    const int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };

    // Add vectors in parallel.
    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}
*/

// Helper function for using CUDA to add vectors in parallel.
/*
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
*/