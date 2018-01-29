#include "JoinTest.h"
#include "gpuKernel.h"

int JoinTest::init(vector<STTraj>* dataPtr, Grid * gridIndex)
{
	this->dataPtr = dataPtr;
	this->gridIndex = gridIndex;
	return 0;
}



int JoinTest::defaultTest(float epsilon, float alpha, int setSize1, int setSize2, map<trajPair, float>& result)
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
	
	this->joinExhaustedCPU(epsilon, alpha, this->taskSet1, this->taskSet2, result);
	timer.stop();
	printf("Elapsed Time: %f ms\n", timer.elapse());
	
	//for (map<trajPair, float>::iterator it = result.begin(); it != result.end(); it++) {
	//	printf("Pair:(%zd, %zd) Distance: %f\t", it->first.first, it->first.second, it->second);
	//}
	
	return 0;
}

int JoinTest::joinExhaustedCPU(float epsilon, float alpha, vector<size_t> &join_set_1, vector<size_t> &join_set_2, map<trajPair, float>& result)
{
	// retrieve two sets of traj
	for (size_t i = 0; i < join_set_1.size(); i++) {
		for (size_t j = 0; j < join_set_2.size(); j++) {
			// calculate similarity distance seperately
			// function: alpha, traj1, traj2 -> dist
			float distance = (*this->dataPtr)[join_set_1[i]].HausdorffDistance((*this->dataPtr)[join_set_2[j]], alpha);
			// verify whether distance is less than epsilon
			if (distance <= epsilon)
				result[trajPair(join_set_1[i], join_set_2[j])] = distance;
		}
	}
	// return result
	return 0;
}

int JoinTest::joinExhaustedGPU(float epsilon, float alpha, vector<size_t>& join_set_1, vector<size_t>& join_set_2, map<trajPair, float>& result)
{

	CUDAwarmUp();
	void* gpuAddr = GPUMalloc((size_t)300 * 1024 * 1024);
	void* gpuAddr8byteAligned = GPUMalloc((size_t)300 * 1024 * 1024);
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	MyTimer timer;

	timer.start();
	for (size_t i = 0; i < join_set_1.size(); i += 255) {
		for (size_t j = 0; j < join_set_2.size(); j += 255) {
			map<trajPair, float> partialResult;
			vector<STTraj> trajSetP, trajSetQ;
			for (size_t i = 0; i < ((i + 255) > join_set_1.size() ? (join_set_1.size()) : (i + 255)); i++) {
				trajSetP.push_back(this->dataPtr->at(join_set_1[i]));
			}
			for (size_t i = 0; i < ((i + 255) > join_set_2.size() ? (join_set_2.size()) : (i + 255)); i++) {
				trajSetQ.push_back(this->dataPtr->at(join_set_2[i]));
			}
			calculateDistanceGPU(trajSetP, trajSetQ, partialResult, gpuAddr, gpuAddr8byteAligned, alpha, epsilon,
				stream, 0);
			// insert new result
			for (map<trajPair, float>::iterator it = partialResult.begin(); it != partialResult.end(); it++) {
				result.insert(*it);
			}
		}
	}

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
