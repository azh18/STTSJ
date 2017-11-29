#include "Cell.h"



Cell::Cell()
{
}

int Cell::isOverlap(MBR mbr)
{
	return bound.isOverlap(mbr);
}

int Cell::insertTraj(STTraj traj)
{
	trajList.insert(traj.trajID);
	return 0;
}

int Cell::getTrajNum()
{
	return trajList.size();
}


Cell::~Cell()
{
}
