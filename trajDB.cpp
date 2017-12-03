#include "trajDB.h"
#include<fstream>
#include<sstream>
#include "gpuKernel.h"
using std::stringstream;
using std::ifstream;
using std::istreambuf_iterator;
using std::string;

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
	// test.defaultTest(epsilon, alpha, setSize1, setSize2);

	// test exhausted join on GPU
	JoinTest testExhaustGPU;
	testExhaustGPU.init(&this->data, &this->gridIndex);
	map<trajPair, double> testResult;
	vector<int> joinsetP, joinsetQ;
	for (int i = 0; i < setSize1; i++)
		joinsetP.push_back(i);
	for (int i = 0; i < setSize2; i++)
		joinsetQ.push_back(i);
	testExhaustGPU.joinExhaustedGPU(epsilon, alpha, joinsetP, joinsetQ, testResult);
	// map<trajPair, double> tt;
	// calculateDistanceGPU((this->data), (this->data), tt);
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
