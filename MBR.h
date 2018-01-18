#pragma once
#include "GlobalHeader.h"
#include "STPoint.h"
class MBR
{
public:
	float lat1, lat2, lon1, lon2;

	MBR();
	MBR(float lat1, float lat2, float lon1, float lon2);
	int isOverlap(MBR mbr2); // return 0:no overlap; 1:mbr2 contain this 2:this contain mbr 3:no contain, only overlap
	int containPoint(STPoint &p);
	~MBR();
};
 
