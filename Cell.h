#pragma once
#include "GlobalHeader.h"
#include "MBR.h"
#include "STTraj.h"

using std::vector;
using std::set;
using std::map;
class Cell
{
public:
	MBR bound;
	int cellid;
	set<size_t> trajList;// pre for expand
	map<size_t, int> trajPointNumList; // map: traid->num
	 
	Cell();
	Cell(int id);
	int isOverlap(MBR mbr);
	size_t insertTraj(STTraj traj);
	size_t getTrajNum();
	~Cell();
};

