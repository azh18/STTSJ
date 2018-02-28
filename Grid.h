#pragma once
#include "GlobalHeader.h"
#include "STPoint.h"
#include "STTraj.h"
#include "Cell.h"
#include <string>


using std::vector;
using std::string;

class Grid
{
public:
	vector<Cell> cells;
	MBR allSpaceBound;
	int gridLatScale, gridLonScale;
	float resl_lat, resl_lon;
	 
	Grid();
	// unittest
	int testAllMethods();
	
	// functions
	int initial(const MBR &overallBound, float resl_lat, float resl_lon);
	int addTrajIntoGrid(const STTraj &traj);
	int getCellIDFromCoordinate(float lat, float lon);
	float getSurroundCellID(const STPoint &p, int probeIter, int cellid, vector<int> &cells);
	int getTrajsOverlappedCell(int cellid, vector<size_t> &trajs);


	//debug
	int outputCellTrajList(string trajListFileName);

	~Grid();
};

