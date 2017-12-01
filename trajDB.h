#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "gpuKernel.h"
#include "Cell.h"
#include "Grid.h"

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

	trajDB();
	int loadDictFromFile(string fileName);
	int loadTrajFromFile(string fileName);
	int cleanOutsideData();
	int buildGridIndex(double resl_lat, double resl_lon);
	int getAllPointNum();

	MBR getMBRofAllData();

	~trajDB();
};

