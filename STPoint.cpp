#include "STPoint.h"

using std::set;
using std::vector;
extern FILE* pfile;

inline bool wordSearch(int wordIdx, vector<int> &wordPool) {
	for (vector<int>::iterator it = wordPool.begin(); it != wordPool.end(); it++) {
		if ((*it) == wordIdx)
			return true;
	}
	return false;
}

STPoint::STPoint(float lat, float lon, std::vector<int> keywords)
{
	this->lat = lat;
	this->lon = lon;
	this->keywords = keywords;
}

STPoint::STPoint(float lat, float lon)
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

float STPoint::STdistance(STPoint & p, float alpha)
{
	// return (this->Sdistance(p)*alpha);
	// fprintf(pfile, "%f\n", this->Tdistance(p));
	return (this->Sdistance(p)*alpha + this->Tdistance(p)*(1 - alpha));
}

float STPoint::Sdistance(STPoint & p)
{
	float lat_dis = this->lat - p.lat;
	float lon_dis = this->lon - p.lon;
	return sqrt(lat_dis*lat_dis + lon_dis*lon_dis) / (float)MAX_DIST;
}

float STPoint::Tdistance(STPoint & p)
{
	// float jaccard;
	//set<int> words;
	//for (size_t i = 0; i < this->keywords.size(); i++) {
	//	words.insert(this->keywords[i]);
	//}
	//for (size_t i = 0; i < p.keywords.size(); i++) {
	//	words.insert(p.keywords[i]);
	//}
	//int intersectSize = 0;
	//for (set<int>::iterator it = words.begin(); it != words.end(); it++) {
	//	if (wordSearch(*it, this->keywords) && wordSearch(*it, p.keywords))
	//		intersectSize++;
	//}
	//jaccard = (float)intersectSize / words.size();
	//return (1.0 - jaccard);

	int tempWords[MAX_KEYWORD_NUM]; uint32_t intersect_size = 0, union_size = 0;
	for (uint32_t idxW = 0; idxW < this->keywords.size(); idxW++) {
		tempWords[idxW] = this->keywords[idxW];
		union_size++;
	}
	for (uint32_t idxW = 0; idxW < p.keywords.size(); idxW++) {
		bool haveSame = false;
		for (uint32_t idxW1 = 0; idxW1 < this->keywords.size(); idxW1++) {
			if (tempWords[idxW1] == p.keywords[idxW]) {
				// intersect_size++;
				haveSame = true;
				break;
			}
		}
		if (haveSame)
			intersect_size++;
		else
			union_size++;
	}
	if (union_size == 0)
		return 0;
	else // do not divide zero!
		return (float)1.0 - (float)intersect_size / union_size;
}



STPoint::~STPoint()
{
}
