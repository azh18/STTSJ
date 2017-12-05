#include "Grid.h"

Grid::Grid()
{
}

int Grid::initial(const MBR & overallBound, double resl_lat, double resl_lon)
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


int Grid::addTrajIntoGrid(const STTraj & traj)
{
	for (size_t i = 0; i < traj.points.size(); i++) {
		int lat_idx = int((traj.points[i].lat - this->allSpaceBound.lat1) / this->resl_lat);
		int lon_idx = int((traj.points[i].lon - this->allSpaceBound.lon1) / this->resl_lon);
		this->cells[lat_idx*this->gridLonScale + lon_idx].insertTraj(traj);
	}
	return 0;
}

// wait for unit test
int Grid::getCellIDFromCoordinate(double lat, double lon)
{
	int lat_idx = (int)((lat - this->allSpaceBound.lat1) / this->resl_lat);
	int lon_idx = (int)((lon - this->allSpaceBound.lon1) / this->resl_lon);
	int cellid = lat_idx * this->gridLonScale + lon_idx;
	return cellid;
}

// get surround cells which are within probeIter cell from cellid
// return the maximum distance of probeIter
// wait for unit test
// if edge is over the overall MBR, this edge will be ignored. (guarantee LB is growing)
double Grid::getSurroundCellID(STPoint &p, int probeIter, int cellid, vector<int> &cells)
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
		cells.push_back(smallestLatIdx * this->gridLonScale + smallestLonIdx);
	}
	// right-up corner
	if (smallestLatIdx >= 0 && largestLonIdx < this->gridLonScale) {
		cells.push_back(smallestLatIdx * this->gridLonScale + smallestLonIdx);
	}
	// right-down corner
	if (largestLatIdx <this->gridLatScale && largestLonIdx < this->gridLonScale) {
		cells.push_back(smallestLatIdx * this->gridLonScale + smallestLonIdx);
	}
	// put edge and compute lowerbound of spatial distance
	double spatialDistanceMin = 9999;
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

int Grid::getTrajsOverlappedCell(int cellid, vector<size_t>& trajs)
{
	for (set<size_t>::iterator it = this->cells[cellid].trajList.begin();
		it != this->cells[cellid].trajList.end(); it++) {
		trajs.push_back(*it);
	}
	return 0;
}


Grid::~Grid()
{
}
