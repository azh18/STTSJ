#include "STPoint.h"

using std::set;
using std::vector;

inline bool wordSearch(int wordIdx, vector<int> &wordPool) {
	for (vector<int>::iterator it = wordPool.begin(); it != wordPool.end(); it++) {
		if ((*it) == wordIdx)
			return true;
	}
	return false;
}

STPoint::STPoint(double lat, double lon, std::vector<int> keywords)
{
	this->lat = lat;
	this->lon = lon;
	this->keywords = keywords;
}

STPoint::STPoint(double lat, double lon)
{
	this->lat = lat;
	this->lon = lon;
}

STPoint::STPoint()
{
	this->lat = 0;
	this->lon = 0;
}

size_t STPoint::getKeywordSize()
{
	return this->keywords.size();
}

double STPoint::STdistance(STPoint & p, double alpha)
{

	return (this->Sdistance(p)*alpha + this->Tdistance(p)*(1 - alpha));
}

double STPoint::Sdistance(STPoint & p)
{
	double lat_dis = this->lat - p.lat;
	double lon_dis = this->lon - p.lon;
	return sqrt(lat_dis*lat_dis + lon_dis*lon_dis) / MAX_DIST;
}

double STPoint::Tdistance(STPoint & p)
{
	double jaccard;
	set<int> words;
	for (size_t i = 0; i < this->keywords.size(); i++) {
		words.insert(this->keywords[i]);
	}
	for (size_t i = 0; i < p.keywords.size(); i++) {
		words.insert(p.keywords[i]);
	}
	int intersectSize = 0;
	for (set<int>::iterator it = words.begin(); it != words.end(); it++) {
		if (wordSearch(*it, this->keywords) && wordSearch(*it, p.keywords))
			intersectSize++;
	}
	jaccard = (double)intersectSize / words.size();
	return (1.0 - jaccard);
}



STPoint::~STPoint()
{
}
