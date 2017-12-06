#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "gpuKernel.h"
#include "Cell.h"
#include "Grid.h"
#include "JoinTest.h"

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
	int cleanOutsideData();
	int buildGridIndex(double resl_lat, double resl_lon);
	size_t getAllPointNum();
	int runDefaultTest(double epsilon, double alpha, int setSize1, int setSize2);
	int getDatasetInformation();
	double similarityGridProber(STPoint &p, set<size_t> &Pset, int probeIter, double alpha, double epsilon,
		set<size_t> &candTrajs, set<size_t> &filteredTrajs,
		vector<map<size_t, bool>> &probedTable, map<size_t, int> &probedTimeArray, int pi);
	int similarityGridFilter(STTraj &t, set<size_t> &Pset, 
		double alpha, double epsilon, 
		vector<size_t> &candTraj);

	MBR getMBRofAllData();
	MBR getMBRofAllData(MBR &mbr);

	//unit test
	int testAllFunctions();



	~trajDB();
};

