#pragma once
#include "GlobalHeader.h"
class STPoint
{
public:
	double lat;
	double lon;
	std::vector<int> keywords;
	STPoint(double lat, double lon, std::vector<int> keywords);
	STPoint();
	int getKeywordSize();
	~STPoint();
};

