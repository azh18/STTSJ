#pragma once
#include "GlobalHeader.h"
class MBR
{
public:
	double xmin, xmax, ymin, ymax;

	MBR();
	MBR(double xmin, double xmax, double ymin, double ymax);
	int isOverlap(MBR mbr2); // return 0:no overlap; 1:mbr2 contain this 2:this contain mbr 3:no contain, only overlap
	~MBR();
};

