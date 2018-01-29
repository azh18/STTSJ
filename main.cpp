#include<iostream>
#include "trajDB.h"

using std::cout;
using std::endl;
FILE* pfile; //write info out

int main() {
	cout << "hello world!" << endl;
	pfile = fopen("outputInfo.txt", "w+");
	STTraj::unitTestForDist();
	trajDB *db = new trajDB();
	db->loadDictFromFile("CPPDict.data");
	db->loadTrajFromFile("trajFile_Twitter.data");
	//int pointNum = db->getAllPointNum();
	//std::cout << pointNum << std::endl;
	db->getMBRofAllData(MBR(39.0, 42.0, -76.0, -72.0));
	db->cleanOutsideData(128);
	db->buildGridIndex((float)0.1, (float)0.1);
	db->buildBloomFilter((float)BLOOM_FILTER_ERROR);
	db->test.init(&db->data, &db->gridIndex);
	db->testAllFunctions((float)0.5, (float)0.2);
	// getchar();
	fclose(pfile);

	// db->getDatasetInformation();
	// db->runDefaultTest((float)0.55,(float)0.5, 128, 128);
	//pointNum = db->getAllPointNum();
	// std::cout << pointNum << std::endl;
	getchar();
	getchar();
	return 0;
}