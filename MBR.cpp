#include "MBR.h"



MBR::MBR()
{
}

MBR::MBR(double xmin, double xmax, double ymin, double ymax)
{
	this->xmin = xmin;
	this->xmax = xmax;
	this->ymin = ymin;
	this->ymax = ymax;
}

int MBR::isOverlap(MBR mbr2)
{
	int overlap;
	if (this->xmax<mbr2.xmin || this->xmin>mbr2.xmax || this->ymin > mbr2.ymax || this->ymax < mbr2.ymin)
		overlap = 0; // not overlap
	else {
		if (this->xmin >= mbr2.xmin && this->xmax <= mbr2.xmax && this->ymin >= mbr2.ymin && this->ymax <= mbr2.ymax)
			overlap = 1; // mbr2 contain this
		else if (this->xmin <= mbr2.xmin && this->xmax >= mbr2.xmax && this->ymin <= mbr2.ymin && this->ymax >= mbr2.ymax)
			overlap = 2; //this contain mbr2
		else
			overlap = 3; // only overlap, no contain
	}
	return overlap;
}


MBR::~MBR()
{
}
