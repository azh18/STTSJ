#include "STPoint.h"

STPoint::STPoint(double x, double y, std::vector<int> keywords)
{
	this->x = x;
	this->y = y;
	this->keywords = keywords;
}

STPoint::STPoint()
{
	this->x = 0;
	this->y = 0;
}

int STPoint::getKeywordSize()
{
	return this->keywords.size();
}


STPoint::~STPoint()
{
}
