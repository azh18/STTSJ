#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "STPoint.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using std::vector;
using std::map;

void* GPUMalloc(size_t byteNum);
int calculateDistanceGPU(const vector<STTraj> &trajSetP,
	const vector<STTraj> &trajSetQ,
	map<trajPair, double> &result,
	void* baseGPUAddr, cudaStream_t &stream);