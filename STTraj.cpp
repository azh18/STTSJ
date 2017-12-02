#include "STTraj.h"
using std::set;
using std::map;

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

size_t STTraj::getLength()
{
	return this->points.size();
}

int STTraj::isOverlap(MBR mbr)
{
	return 0;
}

// unit test: test whether simi calculate is accurate
/*
t1: (1,2)-(3,4),10,11,12
t2: (5,6)-(8-9),11,12,13
distS(p1,p2) = 4sqrt(2), 2sqrt(2), 7sqrt(2), 5sqrt(2)  /6.0
distS(p2,p1) = 4sqrt(2), 2sqrt(2), 7sqrt(2), 5sqrt(2)  /6.0
distT(p1,p2) = distT(p2,p1) = 1-2/4 = 0.5
Hausdorff = 5sqrt(2)/6.0 + 0.5 / 2 = 0.83

*/
int STTraj::unitTestForDist() {
	STTraj t1, t2;
	vector<int> keyword1, keyword2;
	keyword1 = { 10,11,12 };
	keyword2 = { 11,12,13 };
	t1.addPoints(STPoint(1, 2, keyword1));
	t1.addPoints(STPoint(3, 4, keyword1));
	t2.addPoints(STPoint(5, 6, keyword2));
	t2.addPoints(STPoint(8, 9, keyword2));
	vector<STTraj> t1s, t2s;
	printf("Unit Test For Hausdorff Distance: %f\n", t1.HausdorffDistance(t2, 0.5));
	return 0;
}

double STTraj::HausdorffDistance(STTraj & traj, double alpha)
{
	double hausd1, hausd2;
	// compute this to traj
	double maxD = 0;
	for (int i = 0; i < this->points.size(); i++) {
		double minD = 999990;
		for (int j = 0; j < traj.points.size(); j++) {
			double distance = this->points[i].STdistance(traj.points[j], alpha);
			if (distance < minD)
				minD = distance;
		}
		if (minD > maxD)
			maxD = minD;
	}
	hausd1 = maxD;
	// compute traj to this
	maxD = 0;
	for (int i = 0; i < traj.points.size(); i++) {
		double minD = 999990;
		for (int j = 0; j < this->points.size(); j++) {
			double distance = traj.points[i].STdistance(this->points[j], alpha);
			if (distance < minD)
				minD = distance;
		}
		if (minD > maxD)
			maxD = minD;
	}
	hausd2 = maxD;
	return (hausd1 > hausd2 ? hausd1 : hausd2);
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
