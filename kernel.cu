
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "gpuKernel.h"
#include <assert.h>
#include <stdio.h>
#include "GlobalHeader.h"


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

/*
GPU function definition.
All functions of GPU are defined here.

*/
__global__ void addKernel(int *c, const int *a, const int *b)
{
	int i = threadIdx.x;
	c[i] = a[i] + b[i];
}

__device__ inline double SDistance(double lat1, double lat2, double lon1, double lon2) {
	double latdelta = lat1 - lat2;
	double londelta = lon1 - lon2;
	return sqrt(latdelta*latdelta + londelta*londelta);
}

__device__ inline double TDistance(int* word1, int* word2, uint16_t wordNum1, uint16_t wordNum2) {
	return;
}

__global__ void computeHausdorffDistanceByGPU(Latlon* latlonP, int* textP, uint32_t* textIdxPArray,
	Latlon* latlonQ, int* textQ, uint32_t* textIdxQArray, int datasizeP, int datasizeQ, 
	uint16_t *wordNumP, uint16_t *wordNumQ, GPUHausInfoTable* taskInfo) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	GPUHashInfoTable task = taskInfo[bID];
	uint32_t pointIdxP = task.latlonIdxP;
	uint32_t pointNumP = task.pointNumP;
	uint32_t pointIdxQ = task.latlonIdxQ;
	uint32_t pointNumQ = task.pointNumQ;
	// for each thread, get addr and compute distance
	double latP = latlonP[pointIdxP + tID / datasizeQ].lat;
	double lonP = latlonP[pointIdxP + tID / datasizeQ].lon;
	uint32_t textStartIdxP = textIdxPArray[pointIdxP + tID / datasizeQ];
	uint16_t keywordNumP = wordNumP[pointIdxP + tID / datasizeQ];
	int keywordP[MAX_KEYWORD_NUM], keywordQ[MAX_KEYWORD_NUM];
	double latQ = latlonQ[pointIdxQ + tID % datasizeQ].lat;
	double lonQ = latlonQ[pointIdxQ + tID % datasizeQ].lon;
	uint32_t textStartIdxQ = textIdxQArray[pointIdxQ + tID % datasizeQ];
	uint16_t keywordNumQ = wordNumQ[pointIdxQ + tID % datasizeQ];
	__shared__ double distancePt1[THREAD_NUM], distancePt2[THREAD_NUM];



	// 
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
	CUDA_CALL(cudaMalloc(&addr, byteNum));
	return addr;
}


int calculateDistanceGPU(const vector<STTraj> &trajSetP,
	const vector<STTraj> &trajSetQ,
	map<trajPair, double> &result,
	void* baseGPUAddr, cudaStream_t &stream) {
	// Latlon *latlonDataPCPU, *latlonDataQCPU; // latlon array
	// int *textDataPCPU, *textDataQCPU; 
	vector<int> textDataPCPU, textDataQCPU; // keyword array
	vector<uint32_t> textIdxPCPU, textIdxQCPU; // keyword idx for each point (to locate the where and how many keywords for a point)
	vector<uint16_t> numWordPCPU, numWordQCPU; // keyword num in each point
	size_t dataSizeP = trajSetP.size(), dataSizeQ = trajSetQ.size();
	// 所有点线性排列
	vector<Latlon> latlonDataPCPU, latlonDataQCPU;
	GPUHausInfoTable *hausTaskInfoCPU = (GPUHashInfoTable*)malloc(sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ);

	//Latlon *latlon = latlonDataPCPU;
	//size_t keywordNum[1000];
	uint32_t pointCnt = 0, textIdxCnt = 0, textCnt = 0;
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
			numWordPCPU.push_back((uint16_t)trajSetP[i].points[j].keywords.size());
			pointCnt++;
			textIdxPCPU.push_back(textCnt);
			for (size_t k = 0; k < trajSetP[i].points[j].keywords.size(); k++) {
				textDataPCPU.push_back(trajSetP[i].points[j].keywords[k]);
				textCnt++;
			}
		}
	}
	void* latlonDataPGPU, *latlonDataQGPU, *textDataPGPU, *textDataQGPU, *textIdxPGPU, *textIdxQGPU, *numWordPGPU, *numWordQGPU;
	void *pNow = baseGPUAddr;
	// Copy data of P to GPU
	CUDA_CALL(cudaMemcpyAsync(pNow, &latlonDataPCPU[0], sizeof(Latlon)*latlonDataPCPU.size(), cudaMemcpyHostToDevice, stream));
	latlonDataPGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(Latlon)*latlonDataPCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &textDataPCPU[0], sizeof(int)*textDataPCPU.size(), cudaMemcpyHostToDevice, stream));
	textDataPGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(int)*textDataPCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &textIdxPCPU[0], sizeof(uint32_t)*textIdxPCPU.size(), cudaMemcpyHostToDevice, stream));
	textIdxPGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint32_t)*textIdxPCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &numWordPCPU[0], sizeof(uint16_t)*numWordPCPU.size(), cudaMemcpyHostToDevice, stream));
	numWordPGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint16_t)*numWordPCPU.size());

	// process set Q
	pointCnt = 0; textCnt = 0;
	for (size_t i = 0; i < trajSetQ.size(); i++) {
		//update table for Q[][i]
		for (size_t j = i; j < i + dataSizeP*dataSizeQ; j+=dataSizeP) {
			hausTaskInfoCPU[j].pointNumQ = (uint32_t)trajSetQ[i].points.size();
			hausTaskInfoCPU[j].latlonIdxQ = pointCnt;
		}
		//insert data and update idx
		for (size_t j = 0; j < trajSetQ[i].points.size(); j++) {
			Latlon p;
			p.lat = trajSetQ[i].points[j].lat;
			p.lon = trajSetQ[i].points[j].lon;
			latlonDataQCPU.push_back(p);
			numWordQCPU.push_back((uint16_t)trajSetQ[i].points[j].keywords.size());
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
	CUDA_CALL(cudaMemcpyAsync(pNow, &latlonDataQCPU[0], sizeof(Latlon)*latlonDataQCPU.size(), cudaMemcpyHostToDevice, stream));
	latlonDataQGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(Latlon)*latlonDataQCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &textDataQCPU[0], sizeof(int)*textDataQCPU.size(), cudaMemcpyHostToDevice, stream));
	textDataQGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(int)*textDataQCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &textIdxQCPU[0], sizeof(uint32_t)*textIdxQCPU.size(), cudaMemcpyHostToDevice, stream));
	textIdxQGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint32_t)*textIdxQCPU.size());
	CUDA_CALL(cudaMemcpyAsync(pNow, &numWordQCPU[0], sizeof(uint16_t)*numWordQCPU.size(), cudaMemcpyHostToDevice, stream));
	numWordQGPU = pNow;
	pNow = (void*)((char*)pNow + sizeof(uint16_t)*numWordQCPU.size());

	void *taskinfoTable;
	CUDA_CALL(cudaMemcpyAsync(pNow, hausTaskInfoCPU, sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ, cudaMemcpyHostToDevice, stream));
	taskinfoTable = pNow;
	pNow = (void*)((char*)pNow + sizeof(GPUHashInfoTable)*dataSizeP*dataSizeQ);
	// data copy finish, invoke kernel
	double *distanceResult = new double[dataSizeP*dataSizeQ];

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