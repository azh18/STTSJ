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
	// 加载数据，内容通过python代码生成，格式可参考GeoText_Generated文件夹内的py源文件所定义的
	db->loadDictFromFile("CPPDict.data"); // 加载每个word以及其wordID
	db->loadTrajFromFile("trajFile.data"); //加载文本轨迹数据
	//int pointNum = db->getAllPointNum();
	//std::cout << pointNum << std::endl;
	db->getMBRofAllData(MBR(39.0, 42.0, -76.0, -72.0));	// 规定合法坐标范围
	db->cleanOutsideData(128); // 删除不在合法范围内的轨迹点，且仅保留128条轨迹
	db->buildGridIndex((float)0.01, (float)0.01); // 建立grid索引
	db->gridIndex.outputCellTrajList("trajList(cells).txt"); // debug用，将cell的trajList输出
	db->buildBloomFilter((float)BLOOM_FILTER_ERROR); //初始化bloom filter
	db->test.init(&db->data, &db->gridIndex);
	db->testAllFunctions((float)0.5, (float)0.2); // 测试主要功能：索引、相似度计算
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