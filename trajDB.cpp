#include "trajDB.h"
#include<fstream>
#include<sstream>
#include "gpuKernel.h"
using std::stringstream;
using std::ifstream;
using std::istreambuf_iterator;
using std::string;
using std::ofstream;

trajDB::trajDB()
{
}

size_t trajDB::loadDictFromFile(string fileName)
{
	ifstream inputStream(fileName, std::ios_base::in);
	string buffer;
	buffer.assign(istreambuf_iterator<char>(inputStream), istreambuf_iterator<char>());
	stringstream ss;
	ss.str(buffer);
	string linestr;
	while (std::getline(ss, linestr)) {
		size_t commaIdx = linestr.find(',');
		int keywordId = std::atoi(&linestr[commaIdx + 1]);
		string keyword = linestr.substr(0, commaIdx);
		wordDict[keywordId] = keyword;
	}
	return wordDict.size();
}

size_t trajDB::loadTrajFromFile(string fileName)
{
	this->fileName = fileName;
	ifstream trajFS(fileName, std::ios_base::in);
	string ss;
	ss.assign(istreambuf_iterator<char>(trajFS), istreambuf_iterator<char>());
	stringstream buffer;
	buffer.str(ss);
	string linestr;
	while (std::getline(buffer, linestr)) {
		vector<string> pointInfo;
		size_t startfi = 0;
		size_t fi = linestr.find(';');
		STTraj traj;
		while (fi!=linestr.npos)
		{
			pointInfo.push_back(linestr.substr(startfi, fi - startfi));
			startfi = fi + 1;
			fi = linestr.find(';', startfi);
		}
		for (size_t i = 0; i < pointInfo.size(); i++) {
			string pointStr = pointInfo[i];
			size_t keywordIdx = pointStr.find('[');
			size_t latIdx = 0;
			size_t lonIdx = pointStr.find(',') + 1;
			size_t timeIdx = pointStr.find(',', lonIdx) + 1;
			double lat = std::atof(pointStr.substr(latIdx, lonIdx - latIdx - 1).c_str());
			double lon = std::atof(pointStr.substr(lonIdx, timeIdx - lonIdx - 1).c_str());
			size_t nextKeyWord = pointStr.find(',', keywordIdx);
			vector<int> keywords;
			while (nextKeyWord != pointStr.npos) {
				keywords.push_back(std::atoi(pointStr.substr(keywordIdx + 1, nextKeyWord - keywordIdx - 1).c_str()));
				keywordIdx = nextKeyWord;
				nextKeyWord = pointStr.find(',', keywordIdx + 1);
			}
			STPoint pt(lat, lon, keywords);
			traj.addPoints(pt);
		}
		traj.trajID = this->data.size();
		this->data.push_back(traj);
	}
	return this->data.size();
}

int trajDB::cleanOutsideData()
{
	for (size_t i = 0; i < this->data.size(); i++) {
		for (vector<STPoint>::iterator iter = this->data[i].points.begin(); iter!=this->data[i].points.end();) {
			if (this->allDataMBR.containPoint(*iter))
				iter++;
			else
				iter = this->data[i].points.erase(iter);
		}
	}
	return 0;
}

int trajDB::buildGridIndex(double resl_lat, double resl_lon)
{
	// generate Grid using Grid()
	this->gridIndex.initial(this->allDataMBR, resl_lat, resl_lon);
	// insert each traj into Grid
	for (size_t i = 0; i < this->data.size(); i++) {
		this->gridIndex.addTrajIntoGrid(this->data[i]);
	}
	// return success traj num
	return 0;
}

size_t trajDB::getAllPointNum()
{
	size_t totalNum = 0;
	for (size_t i = 0; i < this->data.size(); i++) {
		totalNum += this->data[i].points.size();
	}
	return totalNum;
}

int trajDB::runDefaultTest(double epsilon, double alpha, int setSize1, int setSize2)
{
	// test exhausted join on CPU
	//test.defaultTest(epsilon, alpha, setSize1, setSize2);

	// test exhausted join on GPU
	JoinTest testExhaustGPU;
	testExhaustGPU.init(&this->data, &this->gridIndex);
	map<trajPair, double> testResult;
	vector<size_t> joinsetP, joinsetQ;
	for (size_t i = 0; i < setSize1; i++)
		joinsetP.push_back(i);
	for (size_t i = 0; i < setSize2; i++)
		joinsetQ.push_back(i);
	testExhaustGPU.joinExhaustedGPU(epsilon, alpha, joinsetP, joinsetQ, testResult);
	// map<trajPair, double> tt;
	// calculateDistanceGPU((this->data), (this->data), tt);
	return 0;
}

int trajDB::getDatasetInformation()
{
	int totalWordNum=0;
	map<int, int> wordNumDistribution;
	int totalPointNum=0;
	map<int, int> pointNumDistribution;
	for (size_t i = 0; i < this->data.size(); i++) {
		totalPointNum += (int)this->data[i].points.size();
		if (pointNumDistribution.find((int)this->data[i].points.size()) == pointNumDistribution.end()) {
			pointNumDistribution[(int)this->data[i].points.size()] = 1;
		}
		else
			pointNumDistribution[(int)this->data[i].points.size()]++;
		for (size_t j = 0; j < this->data[i].points.size(); j++) {
			totalWordNum += (int)this->data[i].points[j].keywords.size();
			if (wordNumDistribution.find((int)this->data[i].points[j].keywords.size()) == wordNumDistribution.end()) {
				wordNumDistribution[(int)this->data[i].points[j].keywords.size()] = 1;
			}
			else
				wordNumDistribution[(int)this->data[i].points[j].keywords.size()]++;
		}
	}
	ofstream infoFile("datasetInfo.txt", std::ios_base::out);
	for (map<int, int>::iterator it = wordNumDistribution.begin(); it != wordNumDistribution.end(); it++) {
		infoFile << it->first << ",";
	}
	infoFile << "\n";
	for (map<int, int>::iterator it = wordNumDistribution.begin(); it != wordNumDistribution.end(); it++) {
		infoFile << it->second << ",";
	}
	infoFile << "\n";
	for (map<int, int>::iterator it = pointNumDistribution.begin(); it != pointNumDistribution.end(); it++) {
		infoFile << it->first << ",";
	}
	infoFile << "\n";
	for (map<int, int>::iterator it = pointNumDistribution.begin(); it != pointNumDistribution.end(); it++) {
		infoFile << it->second << ",";
	}
	infoFile.close();
	std::cout << this->data.size() << "trajectories, " <<
		totalPointNum << "points and " << totalWordNum << "words.\n";
	getchar();
	return 0;
}

// probe from a point p in P
// return a lowerbound of unseen trajectories
// wait for unit test
double trajDB::similarityGridProber(STPoint &p, int probeIter, double alpha, double epsilon,
	set<size_t>& candTrajs, set<size_t>& filteredTrajs)
{
	// get cellid of current point p
	int cid_of_p = this->gridIndex.getCellIDFromCoordinate(p.lat, p.lon);
	// find cells should be probed according to probeIter
	vector<int> probedCellIds;
	double lowerboundSpatial = this->gridIndex.getSurroundCellID(p, probeIter, cid_of_p, probedCellIds);
	if (lowerboundSpatial < 0.000001 && probeIter >0) {
		// 说明probe已经超出范围，不必再继续
		return -1.0;
	}
	// get a set of trajs that overlap these cells
	set<size_t> probedTrajs;
	for (size_t i = 0; i < probedCellIds.size(); i++) {
		vector<size_t> trajsInThisCell;
		this->gridIndex.getTrajsOverlappedCell(probedCellIds[i], trajsInThisCell);
		for (int j = 0; j < trajsInThisCell.size(); j++) {
			probedTrajs.insert(trajsInThisCell[j]);
		}
	}
	// for each trajs, if not in two sets, check condition and insert it into sets
	for (set<size_t>::iterator it = probedTrajs.begin(); it != probedTrajs.end(); it++) {
		size_t tid = *it;
		if (candTrajs.find(tid) == candTrajs.end() && filteredTrajs.find(tid) == filteredTrajs.end()) {
			// find dT(min)
			double jaccardMax = 0;
			for (size_t i = 0; i < this->data[tid].points.size(); i++) {
				double tempJaccardMax = min(p.keywords.size(), this->data[tid].points[i].keywords.size()) /
					max(p.keywords.size(), this->data[tid].points[i].keywords.size());
				if (tempJaccardMax > jaccardMax)
					jaccardMax = tempJaccardMax;
			}
			double dTmin = 1 - jaccardMax;
			if (alpha * lowerboundSpatial + (1 - alpha)*dTmin > epsilon)
				filteredTrajs.insert(tid);
			else
				candTrajs.insert(tid);
		}
	}

	return lowerboundSpatial;
}

// wait for unit test
int trajDB::similarityGridFilter(STTraj & t, vector<STTraj>& Pset, double alpha, double epsilon, vector<STTraj>& candTraj)
{
	// for each p in t, probe new trajs
	set<size_t> candTrajSet, filteredTrajSet; // two sets of trajs which is candidate or filtered
	double lowerBoundSpatial = 0;
	int probeIter = 0;
	// start filter phase, until all trajs are distributed in cand and filter, or none of traj can be in cand
	while (candTrajSet.size() + filteredTrajSet.size() < Pset.size() && lowerBoundSpatial <= epsilon / alpha) {
		// calculate the lowerboundSpatial now
		lowerBoundSpatial = 9999.9;
		for (size_t i = 0; i < t.points.size(); i++) {
			// for this point, invoke probe to update candTraj and filteredTraj
			double tempLB = this->similarityGridProber(t.points[i], probeIter,
					alpha, epsilon, candTrajSet, filteredTrajSet);
			if (tempLB < 0) {
				// 无效LB
				continue;
			}
			else {
				lowerBoundSpatial = (lowerBoundSpatial < tempLB ? lowerBoundSpatial : tempLB);
			}
		}
		probeIter++;
	}
	// 结束后所有candTrajSet中的都是candidate，其他的都可以filter掉
	for (set<size_t>::iterator it = candTrajSet.begin(); it != candTrajSet.end(); it++) {
		candTraj.push_back(*it);
	}
	return 0;
}

MBR trajDB::getMBRofAllData()
{
	//double minx=200, miny=200, maxx=-200, maxy=-200;
	//for (size_t i = 0; i < this->data.size(); i++) {
	//	for (size_t j = 0; j < this->data[i].points.size(); j++) {
	//		j = 0;
	//	}
	//}
	this->allDataMBR = MBR(39.0, 42.0, -76.0, -72.0);
	return MBR(39.0, 42.0, -76.0, -72.0);
}




trajDB::~trajDB()
{
}
