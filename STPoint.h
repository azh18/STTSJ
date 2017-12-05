#pragma once
#include "GlobalHeader.h"
class STPoint
{
public:
	double lat;
	double lon;
	std::vector<int> keywords;
	STPoint(double lat, double lon, std::vector<int> keywords);
	STPoint(double lat, double lon);
	STPoint();
	size_t getKeywordSize();
	double STdistance(STPoint &p, double alpha); // wait for construct
	double Sdistance(STPoint &p); //wait for construct
	double Tdistance(STPoint &p); // wait for construct
	~STPoint();
};

 