#include "trajDB.h"
#include<fstream>
#include<sstream>
using std::stringstream;
using std::ifstream;
using std::istreambuf_iterator;
using std::string;

trajDB::trajDB()
{
}

int trajDB::loadDictFromFile(string fileName)
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

int trajDB::loadTrajFromFile(string fileName)
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


trajDB::~trajDB()
{
}
