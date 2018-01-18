#pragma once
#include "GlobalHeader.h"
#include "STPoint.h"
#include "MBR.h"
#include "bloom_filter.hpp"

using std::vector;
class STTraj
{
public:
	vector<STPoint> points;
	size_t trajID;
	float calculateSimilarity(STTraj traj2);
	int addPoints(STPoint p);
	STPoint* getStartPointAddr();
	size_t getLength();
	int isOverlap(MBR mbr);
	float HausdorffDistance(STTraj &traj, float alpha);
	static int unitTestForDist();
	// bloom filter for keywords in this trajectory
	bloom_filter bfWords;
	bloom_parameters bfParam;
	int generateBloomFilter(float errorTole);
	// construction function
	STTraj();
	explicit STTraj(int trajID);
	~STTraj();
};
 
