#pragma once
#include "GlobalHeader.h"
#include "MBR.h"
#include "STTraj.h"

using std::vector;
using std::set;
class Cell
{
public:
	MBR bound;
	int cellid;
	set<int> trajList;

	Cell();
	Cell(int id);
	int isOverlap(MBR mbr);
	int insertTraj(STTraj traj);
	int getTrajNum();
	~Cell();
};

