#pragma once

#include "point.h"

using namespace std;
namespace xmv
{
	//获得点p距离直线lineP1 ~ lineP2的距离
	template <typename T> T GetPointLocation(Point<T> & lineP1, Point<T> & lineP2, Point<T> & p)
	{
		return (lineP1 - p) ^ (lineP2 - p);
	}

	//求两条直线的交点
	template <typename T> Point<T> GetIntersection(Point<T> & line1P1, Point<T> & line1P2, Point<T> & line2P1, Point<T> & line2P2)
	{
		//(y1 - y2)x + (x2 - x1)y + (x1y2 - x2y1) = 0
		T a1 = line1P1.y - line1P2.y;
		T b1 = line1P2.x - line1P1.x;
		T c1 = line1P1.x * line1P2.y - line1P2.x * line1P1.y;

		T a2 = line2P1.y - line2P2.y;
		T b2 = line2P2.x - line2P1.x;
		T c2 = line2P1.x * line2P2.y - line2P2.x * line2P1.y;

		T delta = a1 * b2 - a2 * b1;
		T x = (b1 * c2 - b2 * c1) / delta;
		T y = (c1 * a2 - c2 * a1) / delta;

		return Point<T>(x, y);
	}

	//判断点p是否在三角形p1 ~ p2 ~ p3中
	template <typename T> bool InTriangle(Point<T> & p1, Point<T> & p2, Point<T> & p3, Point<T> & p)
	{
		T v1 = GetPointLocation(p1, p2, p);
		T v2 = GetPointLocation(p2, p3, p);
		T v3 = GetPointLocation(p3, p1, p);
		if (v1 >= 0 && v2 >= 0 && v3 >= 0) return true;
		if (v1 <= 0 && v2 <= 0 && v3 <= 0) return true;
		return false;
	}

	//判断点p是否在四边形p1 ~ p2 ~ p3 ~ p4中
	template <typename T> bool InQuadrilateral(Point<T> & p1, Point<T> & p2, Point<T> & p3, Point<T> & p4,
		Point<T> & p)
	{
		T v1 = GetPointLocation(p1, p2, p);
		T v2 = GetPointLocation(p2, p3, p);
		T v3 = GetPointLocation(p3, p4, p);
		T v4 = GetPointLocation(p4, p1, p);
		if (v1 >= 0 && v2 >= 0 && v3 >= 0 && v4 >= 0) return true;
		if (v1 <= 0 && v2 <= 0 && v3 <= 0 && v4 <= 0) return true;
		return false;
	}

	template <typename T> bool InCheckerboard(Point<T> & p1, Point<T> & p2, Point<T> & p3, Point<T> & p4,
		Point<T> & p)
	{
		T v1 = GetPointLocation(p1, p2, p);
		T v2 = GetPointLocation(p2, p3, p);
		T v3 = GetPointLocation(p3, p4, p);
		T v4 = GetPointLocation(p4, p1, p);
		//if (sign<T>(v1) + sign<T>(v2) + sign<T>(v3) + sign<T>(v4) == 0) return 0;
		return sign<T>(v1) * sign<T>(v2) * sign<T>(v3) * sign<T>(v4) >= 0;
	}
}