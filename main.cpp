#include<iostream>
#include "trajDB.h"

using std::cout;
using std::endl;
int main() {
	cout << "hello world!" << endl;
	trajDB *db = new trajDB();
	db->loadDictFromFile("CPPDict.data");
	db->loadTrajFromFile("trajFile.data");
	getchar();
	getchar();
	return 0;
}