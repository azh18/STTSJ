#include "trajDB.h"
#include<fstream>
#include<sstream>
#include "gpuKernel.h"
#define min(x,y) x>y?y:x
#define max(x,y) x>y?x:y

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
	ifstream fp("DeleteWordList.data", std::ios_base::in);
	std::string buffer;
	set<int> deleteWord;
	while (std::getline(fp, buffer)) {
		int wordidx = std::atoi(buffer.c_str());
		deleteWord.insert(wordidx);
	}

	for (size_t i = 0; i < this->data.size(); i++) {
		for (vector<STPoint>::iterator iter = this->data[i].points.begin(); iter!=this->data[i].points.end();) {
			for (vector<int>::iterator itk = iter->keywords.begin(); itk != iter->keywords.end(); ) {
				if (deleteWord.find(*itk) != deleteWord.end()) {
					itk = iter->keywords.erase(itk);
				}
				else
				{
					itk++;
				}
			}
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

int trajDB::buildBloomFilter(double errorTolerance)
{
	for (int i = 0; i < this->data.size(); i++) {
		this->data[i].generateBloomFilter(errorTolerance);
	}
	return 0;
}

bool doubleEqual(double a, double b) {
	if (a - b < 0.000001 || b - a < 0.000001)
		return true;
	else
		return false;
}
bool compareResult(map<trajPair, double> &result1, map<trajPair, double>& result2,
	vector<trajPair> &differentPair) {
	bool isSame = true;
	if (result1.size() != result2.size())
		return false;
	else {
		for (map<trajPair, double>::iterator it = result1.begin(); it != result1.end(); it++) {
			map<trajPair, double>::iterator it2 = result2.find(it->first);
			if (it2 == result2.end())
			{
				isSame = false;
				differentPair.push_back(it->first);
			}
			else {
				if (!doubleEqual(it->second,it2->second)) {
					isSame = false;
					differentPair.push_back(it->first);
				}
			}
		}
	}
	return isSame;
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
	map<trajPair, double> CPUExaustedResult;
	// test.defaultTest(epsilon, alpha, setSize1, setSize2, CPUExaustedResult);
	// getchar();
	// test exhausted join on GPU
	JoinTest testExhaustGPU;
	testExhaustGPU.init(&this->data, &this->gridIndex);
	map<trajPair, double> GPUExaustedResult;
	vector<size_t> joinsetP, joinsetQ;
	for (size_t i = 0; i < setSize1; i++)
		joinsetP.push_back(i);
	for (size_t i = 0; i < setSize2; i++)
		joinsetQ.push_back(i);

	testExhaustGPU.joinExhaustedGPU(epsilon, alpha, joinsetP, joinsetQ, GPUExaustedResult);
	vector<trajPair> differentPairs;
	if (compareResult(CPUExaustedResult, GPUExaustedResult,differentPairs))
		printf("CPU and GPU get the same result.\n");
	else
		printf("The result of CPU and GPU is not the same!!!\n");
	//vector<size_t> joinsetP1, joinsetQ1;
	//joinsetP1.push_back(0);
	//joinsetQ1.push_back(12);
	//testExhaustGPU.joinExhaustedGPU(epsilon, alpha, joinsetP1, joinsetQ1, testResult);


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
/*
	textualLBType
	0: use original Jaccard max
	1: use lowerbound based on bloom filter
*/
double trajDB::similarityGridProber(STPoint &p, set<size_t> &Pset, int probeIter, double alpha, double epsilon,
	set<size_t>& candTrajs, set<size_t>& filteredTrajs,
	vector<map<size_t, bool>> &probedTable, map<size_t, int> &probedTimeArray, int pi,
	int textualLBType)
{
	// get cellid of current point p
	int cid_of_p = this->gridIndex.getCellIDFromCoordinate(p.lat, p.lon);
	// find cells should be probed according to probeIter
	vector<int> probedCellIds;
	double lowerboundSpatial = this->gridIndex.getSurroundCellID(p, probeIter, cid_of_p, probedCellIds);
	lowerboundSpatial /= MAX_DIST; // normalized
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
			// only process tid that fall into given set
			if (Pset.find(trajsInThisCell[j]) != Pset.end() && 
				probedTable[pi][trajsInThisCell[j]] == false &&
				candTrajs.find(trajsInThisCell[j]) == candTrajs.end() && 
				filteredTrajs.find(trajsInThisCell[j]) == filteredTrajs.end())
			{
				// then it need to be probed...
				probedTrajs.insert(trajsInThisCell[j]);
				probedTable[pi][trajsInThisCell[j]] = true;
			}
		}
	}
	// test
	//printf("\nprobed tid:\n");
	//for (set<size_t>::iterator it = probedTrajs.begin(); it != probedTrajs.end(); it++)
	//	printf("%zd,", *it);
	//printf("\n");
	// test
	// for each trajs, if not in two sets, check condition and insert it into sets
	for (set<size_t>::iterator it = probedTrajs.begin(); it != probedTrajs.end(); it++) {
		size_t tid = *it;
		// find dT(min)
		double jaccardMax = 0;
		if (textualLBType == 0) {
			for (size_t i = 0; i < this->data[tid].points.size(); i++) {
				double tempJaccardMax = (double)(min(p.keywords.size(), this->data[tid].points[i].keywords.size())) /
					max(p.keywords.size(), this->data[tid].points[i].keywords.size());
				if (tempJaccardMax > jaccardMax)
					jaccardMax = tempJaccardMax;
			}
		}
		else if (textualLBType == 1) {
			int maximalMatchNum = 0;
			for (int i = 0; i < p.keywords.size(); i++) {
				if (this->data[tid].bfWords.contains(p.keywords[i]))
					maximalMatchNum++;
			}
			// find minimum union size
			int minimalUnionSize = 0;
			for (int i = 0; i < this->data[tid].points.size(); i++) {
				minimalUnionSize = this->data[tid].points[i].keywords.size() > minimalUnionSize ?
					this->data[tid].points[i].keywords.size() : minimalUnionSize;
			}
			minimalUnionSize += p.keywords.size() - maximalMatchNum;
			if (minimalUnionSize == 0) {
				jaccardMax = 1;
			}
			else
				jaccardMax = maximalMatchNum / minimalUnionSize;
		}
		else {
			std::cerr << "have not define this textualLB type!!!" << std::endl;
		}
		double dTmin = 1 - jaccardMax;
		double dMin = alpha * lowerboundSpatial + (1 - alpha)*dTmin;
			
		if (dMin > epsilon)
		{
			// printf("Traj %zd is filtered. spatial LB is:%f, textual LB is %f, overall LB is %f \n", tid,lowerboundSpatial, dTmin, dMin);
			filteredTrajs.insert(tid);
		}
		else {
			probedTimeArray[tid]--;
			// printf("Insert traj %zd into cand. spatial LB is:%f, textual LB is %f, overall LB is %f\n", tid, lowerboundSpatial, dTmin, dMin);
		}

	}

	return lowerboundSpatial;
}

// return a set of trajID as candidate set
// wait for unit test
int trajDB::similarityGridFilter(STTraj & t, set<size_t>& Pset, double alpha, double epsilon, vector<size_t>& candTraj)
{
	// maintain a table to mark whether for point i, we have known mind(i,j)<epsilon
	vector<map<size_t, bool>> oneConstrainSatisfiedTable(t.points.size());
	for (size_t i = 0; i < t.points.size(); i++)
	{
		for (set<size_t>::iterator it=Pset.begin();it!=Pset.end();it++)
			oneConstrainSatisfiedTable[i][(*it)] = false;
	}
	// maintain an array, reduce 1 when traj is probed by a point, when turn to 0 move to cand set
	map<size_t, int> remainPointNeedToSatisfyArray;
	for (set<size_t>::iterator it = Pset.begin(); it != Pset.end(); it++)
		remainPointNeedToSatisfyArray[*it] = (int)t.points.size();

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
			double tempLB = this->similarityGridProber(t.points[i], Pset, probeIter,
				alpha, epsilon, candTrajSet, filteredTrajSet,
				oneConstrainSatisfiedTable, remainPointNeedToSatisfyArray, (int)i, 0); // need to be adapted
			if (tempLB < 0) {
				// 无效LB
				continue;
			}
			else {
				lowerBoundSpatial = (lowerBoundSpatial < tempLB ? lowerBoundSpatial : tempLB);
			}
		}
		//for (set<size_t>::iterator it = Pset.begin(); it != Pset.end(); it++) {
		//	if (remainPointNeedToSatisfyArray[*it] == 0)
		//		candTrajSet.insert(*it);
		//}
		probeIter++;
		// printf("\n an iteration finish.\n");
	}
	// 结束后所有candTrajSet中的都是candidate，其他的都可以filter掉
	//for (set<size_t>::iterator it = candTrajSet.begin(); it != candTrajSet.end(); it++) {
	//	candTraj.push_back(*it);
	//}
	for (map<size_t, int>::iterator it = remainPointNeedToSatisfyArray.begin();
		it != remainPointNeedToSatisfyArray.end(); it++)
		if (it->second == 0)
			candTraj.push_back(it->first);
	return 0;
}

MBR trajDB::getMBRofAllData()
{
	double minx=200, miny=200, maxx=-200, maxy=-200;
	for (size_t i = 0; i < this->data.size(); i++) {
		for (size_t j = 0; j < this->data[i].points.size(); j++) {
			j = 0;
		}
	}
	return MBR(minx, miny, maxx, maxy);
}

MBR trajDB::getMBRofAllData(MBR &mbr)
{
	this->allDataMBR = mbr;
	return mbr;
}

int trajDB::testAllFunctions()
{
	// test similarityGridProber
	//set<size_t> candTrajs, filteredTrajs, Pset;
	//for (size_t i = 0; i < 99; i++)
	//	Pset.insert(i);
	//vector<map<size_t, bool>> oneConstrainSatisfiedTable(1);
	//for (set<size_t>::iterator it = Pset.begin(); it != Pset.end(); it++) {
	//	oneConstrainSatisfiedTable[0][*it] = false;
	//}
	//map<size_t, int> remainPointNeedToSatisfyArray;
	//for (set<size_t>::iterator it = Pset.begin(); it != Pset.end(); it++) {
	//	remainPointNeedToSatisfyArray[*it] = 1;
	//}
	//this->similarityGridProber(STPoint(39.12, -75.95), Pset, 2, 0.8, 0.15, candTrajs, filteredTrajs,
	//	oneConstrainSatisfiedTable, remainPointNeedToSatisfyArray, 0, 1);


	//printf("test SimilarityGridProber:\n");
	//for(size_t i=0;i<candTrajs.size();i++)
	//	printf("%zd,", *(candTrajs.begin()));
	//printf("\n");
	
	//test similarity filter
	// result: lower bound is too loose
	set<size_t> Qset;
	vector<size_t> candTrajSet;
	for (size_t i = 0; i < 1280; i++)
		Qset.insert(i);
	MyTimer timer;
	timer.start();
	for(size_t i=0;i<1280;i++)
		this->similarityGridFilter(this->data[i], Qset, 0.5, 0.55, candTrajSet);
	timer.stop();
	printf("Filter spend %f ms.\n", timer.elapse());
	printf("Total:%zd,Filtered:%zd\n", Qset.size(), candTrajSet.size());
	//printf("Here are all candidate trajectories:");
	//for (size_t i = 0; i < candTrajSet.size(); i++) {
	//	printf("%zd,", candTrajSet[i]);
	//}
	return 0;
}




trajDB::~trajDB()
{
}
