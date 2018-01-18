#pragma once
#include "GlobalHeader.h"
class STPoint
{
public:
	float lat;
	float lon;
	std::vector<int> keywords;
	STPoint(float lat, float lon, std::vector<int> keywords);
	STPoint(float lat, float lon);
	STPoint();
	size_t getKeywordSize();
	float STdistance(STPoint &p, float alpha); // wait for construct
	float Sdistance(STPoint &p); //wait for construct
	float Tdistance(STPoint &p); // wait for construct
	~STPoint();
};

 