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
	double calculateSimilarity(STTraj traj2);
	int addPoints(STPoint p);
	STPoint* getStartPointAddr();
	size_t getLength();
	int isOverlap(MBR mbr);
	double HausdorffDistance(STTraj &traj, double alpha);
	static int STTraj::unitTestForDist();
	// bloom filter for keywords in this trajectory
	bloom_filter bfWords;
	bloom_parameters bfParam;
	int generateBloomFilter(double errorTole);
	// construction function
	STTraj();
	explicit STTraj(int trajID);
	~STTraj();
};
 
