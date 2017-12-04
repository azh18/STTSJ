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
#define MAX_KEYWORD_NUM 200 //每个点包含keyword的最大值
#define MAX_TRAJ_LENGTH 200
#define THREAD_NUM 256

typedef std::pair<size_t, size_t> trajPair;