#pragma once
#include "GlobalHeader.h"
#include "STTraj.h"
#include "gpuKernel.h"
#include "Cell.h"

using std::vector;
using std::map;
using std::string;
class trajDB
{
public:
	map<int, STTraj> data;
	vector<Cell> gridIndex;
	map<int, string> wordDict;
	string fileName;

	trajDB();
	int loadTrajFromFile(string fileName);
	int buildWordDict();

	~trajDB();
};

