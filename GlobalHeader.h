#pragma once
#include<vector>
#include<set>
#include<map>
#include<iostream>
#include<string>
#include<math.h>
#ifdef WIN32
#include "WinTimer.h"
#else
#include <sys/time.h>
#endif
#ifdef WIN32
#else

 
class MyTimer
{
public:
	MyTimer() {
	};
	double iStart;
	double iEnd;

	double cpuSecond() {
		struct timeval tp;
		gettimeofday(&tp, NULL);
		return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
	}

	inline void start()
	{
		iStart = cpuSecond();
	}
	inline void stop()
	{
		iEnd = cpuSecond();
	}
	inline float elapse()
	{
		return iEnd - iStart;
	}
};
#endif

#define MAX_DIST 6.0 //空间距离的归一化参数

typedef std::pair<size_t, size_t> trajPair;