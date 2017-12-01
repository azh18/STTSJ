#include "STPoint.h"

STPoint::STPoint(double lat, double lon, std::vector<int> keywords)
{
	this->lat = lat;
	this->lon = lon;
	this->keywords = keywords;
}

STPoint::STPoint()
{
	this->lat = 0;
	this->lon = 0;
}

int STPoint::getKeywordSize()
{
	return this->keywords.size();
}


STPoint::~STPoint()
{
}
