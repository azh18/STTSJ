#pragma once
#include "GlobalHeader.h"
#include "STPoint.h"
#include "MBR.h"
using std::vector;
class STTraj
{
public:
	vector<STPoint> points;
	int trajID;
	double calculateSimilarity(STTraj traj2);
	int addPoints(STPoint p);
	STPoint* getStartPointAddr();
	int getLength();
	int isOverlap(MBR mbr);
	STTraj();
	STTraj(int trajID);
	~STTraj();
};

