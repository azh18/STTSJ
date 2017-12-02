#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "STPoint.h"

using std::vector;
using std::map;

int calculateDistanceGPU(vector<STTraj> trajSetP,
	vector<STTraj> trajSetQ,
	map<trajPair, double> &result);