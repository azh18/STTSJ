#include "MBR.h"



MBR::MBR()
{
}

MBR::MBR(double lat1, double lat2, double lon1, double lon2)
{
	this->lat1 = lat1;
	this->lat2 = lat2;
	this->lon1 = lon1;
	this->lon2 = lon2;
}

int MBR::isOverlap(MBR mbr2)
{
	int overlap;
	if (this->lat2<mbr2.lat1 || this->lat1>mbr2.lat2 || this->lon1 > mbr2.lon2 || this->lon2 < mbr2.lon1)
		overlap = 0; // not overlap
	else {
		if (this->lat1 >= mbr2.lat1 && this->lat2 <= mbr2.lat2 && this->lon1 >= mbr2.lon1 && this->lon2 <= mbr2.lon2)
			overlap = 1; // mbr2 contain this
		else if (this->lat1 <= mbr2.lat1 && this->lat2 >= mbr2.lat2 && this->lon1 <= mbr2.lon1 && this->lon2 >= mbr2.lon2)
			overlap = 2; //this contain mbr2
		else
			overlap = 3; // only overlap, no contain
	}
	return overlap;
}

int MBR::containPoint(STPoint & p)
{
	if ((p.lat < this->lat2) && (p.lat > this->lat1) && (p.lon < this->lon2) && (p.lon > this->lon1))
		return true;
	else
		return false;
	return false;
}


MBR::~MBR()
{
}
