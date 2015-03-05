#pragma once

#include <cmath>
using namespace std;

namespace xmv
{

	template<typename T> inline T Abs(T x) { return x >= 0 ? x : -x; }
	template<> inline long Abs<long>(long x) {return x >= 0 ? x : -x; }
	template<> inline int Abs<int>(int x) {return x >= 0 ? x : -x; }
	template<> inline double Abs<double>(double x) {return x >= 0 ? x : -x; }
	template<> inline float Abs<float>(float x) {return x >= 0 ? x : -x; }
	template<> inline short Abs<short>(short x) {return x >= 0 ? x : -x; }
	template<> inline char Abs<char>(char x) {return x >= 0 ? x : -x; }
	template<> inline long long Abs<long long>(long long x) {return x >= 0 ? x : -x; }

	template<typename T> inline T Max(T x, T y) { return x > y ? x : y; }
	template<typename T> inline T Min(T x, T y) { return x < y ? x : y; }

	template<typename T> int Sign(T v) { if (v <= epsilon) return -1;  else if (v >= epsilon) return 1; else return 0;}

	template<typename T> T Sqrt(T x) { return sqrt(x); }

	template<typename T, typename T2> T Pow(T x, T2 y) { return pow(x, y); }

}