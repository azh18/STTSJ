#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "gpuKernel.h"
#include "Cell.h"
#include "Grid.h"
#include "JoinTest.h"
#include "MBR.h"

using std::vector;
using std::map;
using std::string;
class trajDB
{
public:
	vector<STTraj> data;
	Grid gridIndex;
	map<int, string> wordDict;
	string fileName;
	MBR allDataMBR;
	JoinTest test;

	trajDB();
	size_t loadDictFromFile(string fileName);
	size_t loadTrajFromFile(string fileName);
	int cleanOutsideData(int maxDataSize);
	int buildGridIndex(float resl_lat, float resl_lon);
	int buildBloomFilter(float errorTolerance);
	size_t getAllPointNum();
	int runDefaultTest(float epsilon, float alpha, int setSize1, int setSize2);
	int getDatasetInformation();
	float similarityGridProber(STPoint &p, set<size_t> &Pset, int probeIter, float alpha, float epsilon,
		set<size_t> &candTrajs, set<size_t> &filteredTrajs,
		vector<map<size_t, bool>> &probedTable, map<size_t, int> &probedTimeArray, int pi,
		int textualLBType, size_t pTrajID);
	int similarityGridFilter(STTraj &t, set<size_t> &Pset, 
		float alpha, float epsilon, 
		vector<size_t> &candTraj);

	MBR getMBRofAllData();
	MBR getMBRofAllData(const MBR &mbr);

	//unit test
	int testAllFunctions(float alpha, float epsilon);



	~trajDB();
};

