#include<iostream>
#include "trajDB.h"

using std::cout;
using std::endl;


int main() {
	cout << "hello world!" << endl;
	STTraj::unitTestForDist();
	trajDB *db = new trajDB();
	db->loadDictFromFile("CPPDict.data");
	db->loadTrajFromFile("trajFile_Twitter.data");
	//int pointNum = db->getAllPointNum();
	//std::cout << pointNum << std::endl;
	db->getMBRofAllData(MBR(39.0, 42.0, -76.0, -72.0));
	db->cleanOutsideData();
	db->buildGridIndex((float)0.1, (float)0.1);
	db->buildBloomFilter((float)BLOOM_FILTER_ERROR);
	// db->testAllFunctions();
	// getchar();

	db->test.init(&db->data, &db->gridIndex);
	// db->getDatasetInformation();
	db->runDefaultTest((float)0.55,(float)0.5, 128, 128);
	//pointNum = db->getAllPointNum();
	// std::cout << pointNum << std::endl;
	getchar();
	getchar();
	return 0;
}