#pragma once
#include "GlobalHeader.h"
#include "STPoint.h"
#include "STTraj.h"
#include "Cell.h"


using std::vector;

class Grid
{
public:
	vector<Cell> cells;
	MBR allSpaceBound;
	int gridLatScale, gridLonScale;
	double resl_lat, resl_lon;
	 
	Grid();
	// unittest
	int testAllMethods();
	
	// functions
	int initial(const MBR &overallBound, double resl_lat, double resl_lon);
	int addTrajIntoGrid(const STTraj &traj);
	int getCellIDFromCoordinate(double lat, double lon);
	double getSurroundCellID(STPoint &p, int probeIter, int cellid, vector<int> &cells);
	int getTrajsOverlappedCell(int cellid, vector<size_t> &trajs);


	~Grid();
};

