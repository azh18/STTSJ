#include "Grid.h"
#include <set>
#define min(x,y) x>y?y:x
#define max(x,y) x>y?x:y

using std::set;

Grid::Grid()
{
}

int Grid::testAllMethods()
{
	printf("Test getCellIDFromCoordinate:\n");
	printf("Input: scan from cell 1 to cell n\n");
	for (float lat = (float)39.05; lat < 42; lat += (float)0.1) {
		for (float lon = (float)-75.95; lon < -72; lon += (float)0.1) {
			printf("cell id=%d\t", this->getCellIDFromCoordinate(lat, lon));
		}
	}
	printf("Test getSurroundCellID:\n");
	printf("Input: p(39.12,-75.95) at cell(1,0), and probe = 2, so surround=");
	printf("(3,0),(3,1),(3,2),(2,2),(1,2),(0,2), and minDist = 0.15");
	vector<int> surroundCells;
	float minDist = this->getSurroundCellID(STPoint((float)39.12, (float)-75.95), 2, this->gridLonScale, surroundCells);
	printf("output cells:");
	for (int i = 0; i < surroundCells.size(); i++) {
		printf("(%d,%d),", surroundCells[i] / this->gridLonScale, surroundCells[i] % this->gridLonScale);
	}
	printf("\nminDist = %f\n", minDist);

	printf("Test getSurroundCellID:\n");
	printf("Input: p(41.98, -71.95) at cell(29,40), and probe = 2, so surround=");
	printf("(29,37),(28,27),(27,37),(26,37),(26,38),(26,39), and minDist = 0.15");
	surroundCells.clear();
	minDist = this->getSurroundCellID(STPoint((float)41.98, (float)-71.95), 2, this->getCellIDFromCoordinate((float)41.98, (float)-71.95), surroundCells);
	printf("output cells:");
	for (int i = 0; i < surroundCells.size(); i++) {
		printf("(%d,%d),", surroundCells[i] / this->gridLonScale, surroundCells[i] % this->gridLonScale);
	}
	printf("\nminDist = %f\n", minDist);

	return 0;
}

/*
初始化grid索引
*/
int Grid::initial(const MBR & overallBound, float resl_lat, float resl_lon)
{
	this->allSpaceBound = overallBound;
	this->resl_lat = resl_lat;
	this->resl_lon = resl_lon;
	int lat_size, lon_size;
	lat_size = int((overallBound.lat2 - overallBound.lat1) / resl_lat) + 1;
	this->gridLatScale = lat_size;
	lon_size = int((overallBound.lon2 - overallBound.lon1) / resl_lon) + 1;
	this->gridLonScale = lon_size;
	for (int i = 0; i < lat_size; i++) {
		for (int j = 0; j < lon_size; j++) {
			Cell c;
			c.bound = MBR(overallBound.lat1 + i*resl_lat, overallBound.lat1 + (i + 1)*resl_lat,
				overallBound.lon1 + j*resl_lon, overallBound.lon1 + (j + 1)*resl_lon);
			c.cellid = i*lon_size + j;
			this->cells.push_back(c);
		}
	}
	return 0;
}

/*
将连续存储的轨迹数据导入grid中
*/
int Grid::addTrajIntoGrid(const STTraj & traj)
{
	for (size_t i = 0; i < traj.points.size(); i++) {
		int lat_idx = int((traj.points[i].lat - this->allSpaceBound.lat1) / this->resl_lat);
		int lon_idx = int((traj.points[i].lon - this->allSpaceBound.lon1) / this->resl_lon);
		this->cells[lat_idx*this->gridLonScale + lon_idx].insertTraj(traj);
	}
	return 0;
}

/*
在row-major coding下，给出坐标经纬度，确定其位于的cell的id
*/
// test pass
int Grid::getCellIDFromCoordinate(float lat, float lon)
{
	int lat_idx = (int)((lat - this->allSpaceBound.lat1) / this->resl_lat);
	int lon_idx = (int)((lon - this->allSpaceBound.lon1) / this->resl_lon);
	int cellid = lat_idx * this->gridLonScale + lon_idx;
	return cellid;
}

// get surround cells which are within probeIter cell from cellid
// return the maximum distance of probeIter
// unit test pass
// if edge is over the overall MBR, this edge will be ignored. (guarantee LB is growing)
float Grid::getSurroundCellID(const STPoint &p, int probeIter, int cellid, vector<int> &cells)
{
	int lat_idx = cellid / this->gridLonScale;
	int lon_idx = cellid % this->gridLonScale;
	// consider edge code
	//int smallestLatIdx = (lat_idx - probeIter) < 0 ? 0 : (lat_idx - probeIter);
	//int largestLatIdx = (lat_idx + probeIter) > (this->gridLatScale - 1) ? (this->gridLatScale - 1) : (lat_idx + probeIter);
	//int smallestLonIdx = (lon_idx - probeIter) < 0 ? 0 : (lon_idx - probeIter);
	//int largestLonIdx = (lon_idx + probeIter) > (this->gridLonScale - 1) ? (this->gridLonScale - 1) : (lon_idx + probeIter);
	//cells.push_back(smallestLatIdx * this->gridLonScale + smallestLonIdx);
	//if (largestLonIdx != smallestLonIdx)
	//	cells.push_back(smallestLatIdx * this->gridLonScale + largestLonIdx);
	//if (largestLatIdx != smallestLonIdx) {
	//	cells.push_back(largestLatIdx * this->gridLonScale + smallestLonIdx);
	//	if (largestLonIdx != smallestLonIdx)
	//		cells.push_back(largestLatIdx * this->gridLonScale + largestLonIdx);
	//}
	//// push inside in //error!!!!!
	//for (int i = smallestLatIdx + 1; i < largestLatIdx; i++) {
	//	for (int j = smallestLonIdx + 1; j < largestLonIdx; j++) {
	//		cells.push_back(i*this->gridLonScale + j);
	//	}
	//}
	//smallestLatIdx = (lat_idx - probeIter + 1) < 0 ? 0 : (lat_idx - probeIter + 1);
	//largestLatIdx = (lat_idx + probeIter - 1) > (this->gridLatScale - 1) ? (this->gridLatScale - 1) : (lat_idx + probeIter - 1);
	//smallestLonIdx = (lon_idx - probeIter + 1) < 0 ? 0 : (lon_idx - probeIter + 1);
	//largestLonIdx = (lon_idx + probeIter - 1) > (this->gridLonScale - 1) ? (this->gridLonScale - 1) : (lon_idx + probeIter - 1);
	//// calculate last iter corner coordinate
	//double smallestLon = this->allSpaceBound.lon1 + smallestLonIdx*this->resl_lon;
	//double smallestLat = this->allSpaceBound.lat1 + smallestLatIdx*this->resl_lat;
	//double largestLon = this->allSpaceBound.lon1 + (largestLonIdx + 1)*this->resl_lon;
	//double largestLat = this->allSpaceBound.lat1 + (largestLatIdx + 1)*this->resl_lat;
	//// 
	//return maxDist;

	// ignore edge code
	int smallestLatIdx = (lat_idx - probeIter);
	int largestLatIdx = (lat_idx + probeIter);
	int smallestLonIdx = (lon_idx - probeIter);
	int largestLonIdx = (lon_idx + probeIter);
	// push corner in
	// left-up corner
	if (smallestLatIdx >= 0 && smallestLonIdx >=0) {
		cells.push_back(smallestLatIdx * this->gridLonScale + smallestLonIdx);
	}
	// left-down corner
	if (largestLatIdx < this->gridLatScale && smallestLonIdx >= 0) {
		cells.push_back(largestLatIdx * this->gridLonScale + smallestLonIdx);
	}
	// right-up corner
	if (smallestLatIdx >= 0 && largestLonIdx < this->gridLonScale) {
		cells.push_back(smallestLatIdx * this->gridLonScale + largestLonIdx);
	}
	// right-down corner
	if (largestLatIdx <this->gridLatScale && largestLonIdx < this->gridLonScale) {
		cells.push_back(largestLatIdx * this->gridLonScale + largestLonIdx);
	}
	// put edge and compute lowerbound of spatial distance
	float spatialDistanceMin = 9999;
	// up edge
	if (smallestLatIdx >= 0) {
		for (int i = (smallestLonIdx >= 0 ? smallestLonIdx + 1 : 0);
			i <= (largestLonIdx < this->gridLonScale ? largestLonIdx - 1 : this->gridLonScale - 1);
			i++) {
			cells.push_back(smallestLatIdx*this->gridLonScale + i);
		}
		spatialDistanceMin = min(spatialDistanceMin, p.lat - (this->allSpaceBound.lat1 + (smallestLatIdx + 1) *this->resl_lat));
	}
	// down edge
	if (largestLatIdx < this->gridLatScale) {
		for (int i = (smallestLonIdx >= 0 ? smallestLonIdx + 1 : 0);
			i <= (largestLonIdx < this->gridLonScale ? largestLonIdx - 1 : this->gridLonScale - 1);
			i++) {
			cells.push_back(largestLatIdx*this->gridLonScale + i);
		}
		spatialDistanceMin = min(spatialDistanceMin, (this->allSpaceBound.lat1 + (largestLatIdx)*this->resl_lat) - p.lat);
	}
	// left edge
	if (smallestLonIdx >= 0) {
		for (int i = (smallestLatIdx >= 0 ? smallestLatIdx + 1 : 0);
			i <= (largestLatIdx < this->gridLatScale ? largestLatIdx - 1 : this->gridLatScale - 1);
			i++) {
			cells.push_back(i*this->gridLonScale + smallestLonIdx);
		}
		spatialDistanceMin = min(spatialDistanceMin, p.lon - (this->allSpaceBound.lon1 + (smallestLonIdx + 1)*this->resl_lon));
	}
	//right edge
	if (largestLonIdx < this->gridLonScale) {
		for (int i = (smallestLatIdx >= 0 ? smallestLatIdx + 1 : 0);
			i <= (largestLatIdx < this->gridLatScale ? largestLatIdx - 1 : this->gridLatScale - 1);
			i++) {
			cells.push_back(i*this->gridLonScale + largestLonIdx);
		}
		spatialDistanceMin = min(spatialDistanceMin,(this->allSpaceBound.lon1 + largestLonIdx * this->resl_lon) - p.lon);
	}


	// calculate last iter corner coordinate
	if (probeIter == 0 || spatialDistanceMin == 9999)
		return 0.0;
	else
		return spatialDistanceMin;


}

/*
查找与某个cell所相交的轨迹的编号，放进trajs中
*/
int Grid::getTrajsOverlappedCell(int cellid, vector<size_t>& trajs)
{
	for (set<size_t>::iterator it = this->cells[cellid].trajList.begin();
		it != this->cells[cellid].trajList.end(); it++) {
		trajs.push_back(*it);
	}
	return 0;
}

/*
debug用：输出所有cell所被交叉的轨迹到文件中
*/
int Grid::outputCellTrajList(string trajListFileName)
{
	FILE* p = fopen(trajListFileName.c_str(), "w+");
	fprintf(p, "%d,%d\n", this->gridLatScale, this->gridLonScale);
	for (int i = 0; i < this->cells.size(); i++) {
		fprintf(p, "%d:", i); // the i-th cell
		for (set<size_t>::iterator it = this->cells[i].trajList.begin(); it != this->cells[i].trajList.end(); it++) {
			fprintf(p, "%d,", *it);
		}
		fprintf(p, "\n");
	}
	return 0;
}

Grid::~Grid()
{
}
