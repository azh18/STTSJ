
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "gpuKernel.h"

#include <stdio.h>

typedef struct GPUHausInfoTable {
	size_t *keywordNumP, *keywordNumQ; // #keywords in each point
	size_t taskNumP, taskNumQ; // #traj in each set
	size_t *pointNumP, *pointNumQ; // #points in each traj
};

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

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


int calculateDistanceGPU(vector<STTraj> trajSetP,
	vector<STTraj> trajSetQ,
	map<trajPair, double> &result) {
	char *dataSetP, *dataSetQ;
	size_t dataSizeP = 0, dataSizeQ = 0;
	// 所有点线性排列
	dataSizeP = calculateDatasize_TrajSet(trajSetP);
	dataSizeQ = calculateDatasize_TrajSet(trajSetQ);
	dataSetP = (char*)malloc(dataSizeP);
	dataSetQ = (char*)malloc(dataSizeQ);
	GPUHausInfoTable hausTaskInfo;
	

	vector<size_t> pointNumPCPU, pointNumQCPU,keywordNumPCPU, keywordNumQCPU;

	char *p = dataSetP, *q = dataSetQ;
	size_t copiedDataSize = 0;
	size_t keywordNum[1000];
	for (vector<STTraj>::iterator it = trajSetP.begin(); it != trajSetP.end(); it++) {
		copiedDataSize = copySTTrajToArray(*it, p, keywordNum);
		for (int i = 0; i < it->points.size(); i++) {
			keywordNumPCPU.push_back(keywordNum[i]);
		}
		pointNumPCPU.push_back(it->points.size());
		p = p + copiedDataSize;
	}
	for (vector<STTraj>::iterator it = trajSetQ.begin(); it != trajSetQ.end(); it++) {
		copiedDataSize = copySTTrajToArray(*it, q, keywordNum);
		for (int i = 0; i < it->points.size(); i++) {
			keywordNumQCPU.push_back(keywordNum[i]);
		}
		pointNumQCPU.push_back(it->points.size());
		q = q + copiedDataSize;
	}

	return 0;




}


__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
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
