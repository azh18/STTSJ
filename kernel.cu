
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
	// size_t textIdxP, textIdxQ; //the first offset of text of the traj 
	uint32_t pointNumP, pointNumQ; // # point in each traj
	// uint32_t idxInTextIdxP, idxInTextIdxQ; //not need, because is same as latlonIdx// idx of first point in textIdx, from this to derive the range of keywords of each point
}GPUHashInfoTable;

typedef struct Latlon {
	double lat;
	double lon;
}Latlon;

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);
void CUDAwarmUp() {
	CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));
	CUDA_CALL(cudaSetDevice(0));
}

static __inline__ __device__ double atomicMin(double* address, double val) {
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address, assumed;
	assumed = *address_as_ull;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val > __longlong_as_double(assumed) ?
			__longlong_as_double(assumed) : val));
	} while (assumed != old);
	return __longlong_as_double(old);
}

/*
GPU function definition.
All functions of GPU are defined here.

*/
__global__ void addKernel(int *c, const int *a, const int *b)
{
	int i = threadIdx.x;
	c[i] = a[i] + b[i];
}

__device__ inline double SDistance(double lat1, double lon1, double lat2, double lon2) {
	double latdelta = lat1 - lat2;
	double londelta = lon1 - lon2;
	return sqrt(latdelta*latdelta + londelta*londelta)/ MAX_DIST;
}

__device__ inline double TDistance(int* word1, int* word2, uint32_t wordNum1, uint32_t wordNum2) {
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
	return 1.0 - (double)intersect_size / union_size;
}

__global__ void computeHausdorffDistanceByGPU(Latlon* latlonP, int* textP, uint32_t* textIdxPArray,
	Latlon* latlonQ, int* textQ, uint32_t* textIdxQArray, int datasizeP, int datasizeQ, 
	uint32_t *wordNumP, uint32_t *wordNumQ, GPUHausInfoTable* taskInfo, double alpha,
	double *HausdorffResult) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	GPUHashInfoTable task = taskInfo[bID];
	uint32_t pointIdxP = task.latlonIdxP;
	uint32_t pointNumP = task.pointNumP;
	uint32_t pointIdxQ = task.latlonIdxQ;
	uint32_t pointNumQ = task.pointNumQ;
	// for each thread, get addr and compute distance
	// if traj length is larger than 16, use for loop to process
	__shared__ double tempDist[THREAD_NUM];
	__shared__ double minimunDist[MAX_TRAJ_LENGTH+2]; // one more is because easy to reduce maximum
	__shared__ double maxDist1;

	if (tID < pointNumP)
		minimunDist[tID] = 9999.0;
	else if (tID < MAX_TRAJ_LENGTH + 2)
		minimunDist[tID] = -1.0;	
	double latP, latQ, lonP, lonQ;
	uint32_t textStartIdxP, textStartIdxQ;
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
				double TDist = TDistance(&textP[textStartIdxP], &textQ[textStartIdxQ], keywordNumP, keywordNumQ);
				double SDist = SDistance(latQ, lonQ, latP, lonP);
				// dist = S*alpha + T*(1-alpha)
				tempDist[tID] = SDist * alpha + TDist * (1 - alpha);
			}
			// compute minimum, loop (may optimize to reduce)
			// only use a warp
			__syncthreads();
			double minDist = 99999;
			if (tID / 16 == 0 && (tID % 16) < (pointNumP - i)) {
				for (int k = 0; k + j < pointNumQ && k < 16; k++) {
					if (tempDist[k + tID * 16] < minDist)
						minDist = tempDist[k + tID * 16];
				}
				
				//if (minimunDist[tID % 16 + i] > minDist)
				//	minimunDist[tID % 16 + i] = minDist;
				atomicMin(&minimunDist[tID % 16 + i], minDist);
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
				double TDist = TDistance(&textQ[textStartIdxQ], &textP[textStartIdxP], keywordNumQ, keywordNumP);
				double SDist = SDistance(latQ, lonQ, latP, lonP);
				// dist = S*alpha + T*(1-alpha)
				tempDist[tID] = SDist * alpha + TDist * (1 - alpha);
			}
			// compute minimum, loop (may optimize to reduce)
			// only use a warp
			__syncthreads();
			double minDist = 9999999;
			if (tID / 16 == 0 && (tID % 16) < (pointNumQ - i)) {
				for (int k = 0; k + j < pointNumP && k < 16; k++) {
					if (tempDist[k + tID * 16] < minDist)
						minDist = tempDist[k + tID * 16];
				}
				
				atomicMin(&minimunDist[tID % 16 + i], minDist);
				//if (minimunDist[tID % 16 + i] > minDist)
				//	minimunDist[tID % 16 + i] = minDist;
			}
			__syncthreads();
		}
	}
	// 归并找最大值
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

/*
CPU function definition.
All functions of CPU are defined here.

*/

size_t calculateDatasize_TrajSet(vector<STTraj> &trajSet) {
	size_t datasize = 0;
	for (vector<STTraj>::iterator it = trajSet.begin(); it != trajSet.end(); it++) {
		for (vector<STPoint>::iterator itp = it->points.begin(); itp != it->points.end(); itp++) {
			datasize += (2 * sizeof(double) + (itp->keywords.size()) * sizeof(int));
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
		memcpy(pStart, &itp->lat, sizeof(double));
		pStart += sizeof(double);
		memcpy(pStart, &itp->lon, sizeof(double));
		pStart += sizeof(double);
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
	map<trajPair, double> &result,
	void* baseGPUAddr, void* baseGPUAddr8byteAligned, double alpha, double epsilon,
	cudaStream_t &stream) {
	// Latlon *latlonDataPCPU, *latlonDataQCPU; // latlon array
	// int *textDataPCPU, *textDataQCPU; 
	vector<int> textDataPCPU, textDataQCPU; // keyword array
	vector<uint32_t> textIdxPCPU, textIdxQCPU; // keyword idx for each point (to locate the where and how many keywords for a point)
	vector<uint32_t> numWordPCPU, numWordQCPU; // keyword num in each point
	size_t dataSizeP = trajSetP.size(), dataSizeQ = trajSetQ.size();
	// 所有点线性排列
	vector<Latlon> latlonDataPCPU, latlonDataQCPU;
	GPUHausInfoTable *hausTaskInfoCPU = (GPUHashInfoTable*)malloc(sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ);

	//MyTimer timer;
	//timer.start();
	//Latlon *latlon = latlonDataPCPU;
	//size_t keywordNum[1000];
	uint32_t pointCnt = 0,  textCnt = 0;
	// process set P
	for (size_t i = 0; i < trajSetP.size(); i++) {
		//update table for P[i][]
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
			for (size_t k = 0; k < trajSetP[i].points[j].keywords.size(); k++) {
				textDataPCPU.push_back(trajSetP[i].points[j].keywords[k]);
				textCnt++;
			}
		}
	}
	void* latlonDataPGPU, *latlonDataQGPU, *textDataPGPU, *textDataQGPU, *textIdxPGPU, *textIdxQGPU, *numWordPGPU, *numWordQGPU;
	void *pNow = baseGPUAddr, *pNow8ByteAligned = baseGPUAddr8byteAligned;
	// Copy data of P to GPU
	CUDA_CALL(cudaMemcpyAsync(pNow8ByteAligned, &latlonDataPCPU[0], sizeof(Latlon)*latlonDataPCPU.size(), cudaMemcpyHostToDevice, stream));
	latlonDataPGPU = pNow8ByteAligned;
	pNow8ByteAligned = (void*)((char*)pNow8ByteAligned + sizeof(Latlon)*latlonDataPCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &textDataPCPU[0], sizeof(int)*textDataPCPU.size(), cudaMemcpyHostToDevice, stream));
	textDataPGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(int)*textDataPCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &textIdxPCPU[0], sizeof(uint32_t)*textIdxPCPU.size(), cudaMemcpyHostToDevice, stream));
	textIdxPGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint32_t)*textIdxPCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &numWordPCPU[0], sizeof(uint32_t)*numWordPCPU.size(), cudaMemcpyHostToDevice, stream));
	numWordPGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint32_t)*numWordPCPU.size());

	// process set Q
	pointCnt = 0; textCnt = 0;
	for (size_t i = 0; i < trajSetQ.size(); i++) {
		//update table for Q[][i]
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
			pointCnt++;
			textIdxQCPU.push_back(textCnt);
			for (size_t k = 0; k < trajSetQ[i].points[j].keywords.size(); k++) {
				textDataQCPU.push_back(trajSetQ[i].points[j].keywords[k]);
				textCnt++;
			}
		}
	}
	// cpu data build finished
	// transfer data to GPU
	// Copy data of P to GPU
	// this order is to guarantee aligned load in GPU

	CUDA_CALL(cudaMemcpyAsync(pNow8ByteAligned, &latlonDataQCPU[0], sizeof(Latlon)*latlonDataQCPU.size(), cudaMemcpyHostToDevice, stream));
	latlonDataQGPU = pNow8ByteAligned;
	pNow8ByteAligned = (void*)((char*)pNow8ByteAligned + sizeof(Latlon)*latlonDataQCPU.size());

	CUDA_CALL(cudaMemcpyAsync(pNow, &textDataQCPU[0], sizeof(int)*textDataQCPU.size(), cudaMemcpyHostToDevice, stream));
	textDataQGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(int)*textDataQCPU.size());

	CUDA_CALL(cudaMemcpyAsync(pNow, &textIdxQCPU[0], sizeof(uint32_t)*textIdxQCPU.size(), cudaMemcpyHostToDevice, stream));
	textIdxQGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint32_t)*textIdxQCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &numWordQCPU[0], sizeof(uint32_t)*numWordQCPU.size(), cudaMemcpyHostToDevice, stream));
	numWordQGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint32_t)*numWordQCPU.size());

	void *taskinfoTable;
	CUDA_CALL(cudaMemcpyAsync(pNow, hausTaskInfoCPU, sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ, cudaMemcpyHostToDevice, stream));
	taskinfoTable = pNow;
	pNow = (void*)((char*)pNow + sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ);
	// data copy finish, invoke kernel
	double *distanceResult,*distanceResultGPU;
	CUDA_CALL(cudaHostAlloc((void**)&distanceResult, dataSizeP*dataSizeQ*sizeof(double), cudaHostAllocMapped));
	CUDA_CALL(cudaHostGetDevicePointer((void**)&distanceResultGPU, distanceResult, 0));
	//cudaDeviceSynchronize();
	//timer.stop();
	//std::cout << timer.elapse() << "ms to copy data" << std::endl;
	computeHausdorffDistanceByGPU << <(uint32_t)dataSizeP*(uint32_t)dataSizeQ, THREAD_NUM, 0, stream >> > ((Latlon*)latlonDataPGPU, 
		(int*)textDataPGPU, (uint32_t*)textIdxPGPU, (Latlon*)latlonDataQGPU, 
		(int*)textDataQGPU, (uint32_t*)textIdxQGPU, 
		(int)dataSizeP, (int)dataSizeQ,
		(uint32_t*)numWordPGPU, (uint32_t*)numWordQGPU, 
		(GPUHashInfoTable*)taskinfoTable, alpha,
		(double*)distanceResultGPU);
	cudaDeviceSynchronize();
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