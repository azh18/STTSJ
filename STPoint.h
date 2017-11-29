#pragma once
#include "GlobalHeader.h"
class STPoint
{
public:
	double x;
	double y;
	std::vector<int> keywords;
	STPoint(double x, double y, std::vector<int> keywords);
	STPoint();
	int getKeywordSize();
	~STPoint();
};

