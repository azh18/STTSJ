#include "STTraj.h"



double STTraj::calculateSimilarity(STTraj traj2)
{
	return 0.0;
}

int STTraj::addPoints(STPoint p)
{
	points.push_back(p);
	return 0;
}

STPoint * STTraj::getStartPointAddr()
{
	if (points.size() > 0)
		return &points[0];
	else
		return NULL;
}

int STTraj::getLength()
{
	return this->points.size();
}

int STTraj::isOverlap(MBR mbr)
{
	return 0;
}

STTraj::STTraj()
{
}

STTraj::STTraj(int trajID)
{
	this->trajID = trajID;
}


STTraj::~STTraj()
{
}
