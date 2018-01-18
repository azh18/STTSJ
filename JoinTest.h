#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "Grid.h"


using std::vector;
using std::map;
class JoinTest
{
public:
	vector<STTraj> *dataPtr;
	Grid *gridIndex;
	vector<size_t> taskSet1, taskSet2;
	int init(vector<STTraj> *dataPtr, Grid *gridIndex);
	int defaultTest(float epsilon, float alpha, int setSize1, int setSize2, map<trajPair, float>& result);
	// join operations...
	int joinExhaustedCPU(float epsilon,
		float alpha,
		vector<size_t> &join_set_1,
		vector<size_t> &join_set_2,
		map<trajPair, float> &result);
	int joinExhaustedGPU(float epsilon,
		float alpha,
		vector<size_t> &join_set_1,
		vector<size_t> &join_set_2,
		map<trajPair, float> &result);

	JoinTest();
	~JoinTest();
};

 