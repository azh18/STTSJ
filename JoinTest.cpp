#include "JoinTest.h"
#include "gpuKernel.h"

int JoinTest::init(vector<STTraj>* dataPtr, Grid * gridIndex)
{
	this->dataPtr = dataPtr;
	this->gridIndex = gridIndex;
	return 0;
}



int JoinTest::defaultTest(double epsilon, double alpha, int setSize1, int setSize2)
{
	printf("starting...\n");

	MyTimer timer;

	timer.start();
	for (size_t i = 0; i < setSize1; i++) {
		this->taskSet1.push_back(i);
	}
	for (size_t i = 0; i < setSize2; i++) {
		this->taskSet2.push_back(i);
	}
	map<trajPair, double> result;
	this->joinExhaustedCPU(epsilon, alpha, this->taskSet1, this->taskSet2, result);
	timer.stop();
	printf("Elapsed Time: %f ms\n", timer.elapse());
	
	for (map<trajPair, double>::iterator it = result.begin(); it != result.end(); it++) {
		printf("Pair:(%zd, %zd) Distance: %f\t", it->first.first, it->first.second, it->second);
	}
	
	return 0;
}

int JoinTest::joinExhaustedCPU(double epsilon, double alpha, vector<size_t> &join_set_1, vector<size_t> &join_set_2, map<trajPair, double>& result)
{
	// retrieve two sets of traj
	for (size_t i = 0; i < join_set_1.size(); i++) {
		for (size_t j = 0; j < join_set_2.size(); j++) {
			// calculate similarity distance seperately
			// function: alpha, traj1, traj2 -> dist
			double distance = (*this->dataPtr)[join_set_1[i]].HausdorffDistance((*this->dataPtr)[join_set_2[j]], alpha);
			// verify whether distance is less than epsilon
			if (distance <= epsilon)
				result[trajPair(join_set_1[i], join_set_2[j])] = distance;
		}
	}
	// return result
	return 0;
}

int JoinTest::joinExhaustedGPU(double epsilon, double alpha, vector<size_t>& join_set_1, vector<size_t>& join_set_2, map<trajPair, double>& result)
{
	vector<STTraj> trajSetP, trajSetQ;
	for (size_t i = 0; i < join_set_1.size(); i++) {
		trajSetP.push_back(this->dataPtr->at(join_set_1[i]));
	}
	for (size_t i = 0; i < join_set_2.size(); i++) {
		trajSetQ.push_back(this->dataPtr->at(join_set_2[i]));
	}
	CUDAwarmUp();
	void* gpuAddr = GPUMalloc((size_t)300 * 1024 * 1024);
	void* gpuAddr8byteAligned = GPUMalloc((size_t)300 * 1024 * 1024);
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	MyTimer timer;

	timer.start();

	calculateDistanceGPU(trajSetP, trajSetQ, result, gpuAddr, gpuAddr8byteAligned, alpha, epsilon,
		stream);
	timer.stop();
	printf("Elapsed Time(GPU Exhausted): %f ms\n", timer.elapse());
	// calculateDistanceGPU(trajSetP, trajSetQ, result);
	return 0;
}

JoinTest::JoinTest()
{
}


JoinTest::~JoinTest()
{
}
