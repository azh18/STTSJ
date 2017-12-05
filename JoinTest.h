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
	int defaultTest(double epsilon, double alpha, int setSize1, int setSize2);
	// join operations...
	int joinExhaustedCPU(double epsilon,
		double alpha,
		vector<size_t> &join_set_1,
		vector<size_t> &join_set_2,
		map<trajPair, double> &result);
	int joinExhaustedGPU(double epsilon,
		double alpha,
		vector<size_t> &join_set_1,
		vector<size_t> &join_set_2,
		map<trajPair, double> &result);

	JoinTest();
	~JoinTest();
};

 