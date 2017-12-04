#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "STPoint.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using std::vector;
using std::map;

void CUDAwarmUp();
void* GPUMalloc(size_t byteNum);
int calculateDistanceGPU(vector<STTraj> &trajSetP,
	vector<STTraj> &trajSetQ,
	map<trajPair, double> &result,
	void* baseGPUAddr, double alpha, cudaStream_t &stream);