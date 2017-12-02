#pragma once
#include "GlobalHeader.h"
#include "STPoint.h"
#include "MBR.h"
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
	STTraj();
	STTraj(int trajID);
	~STTraj();
};
 
