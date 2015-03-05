#pragma once

#include <istream>
#include <cmath>
#include "convert.h"

namespace xmv
{
	//平面中的一个点
	template <typename T> class Point
	{
	public:
		T x;
		T y;

		Point() { }
		Point(T _x, T _y) { x = _x; y = _y; }
		Point(T _x) { x = _x; y = 0; }

		T Modulo() { return data(x * x + y * y); }
		double Angle(){ return atan2(y, x); }
		T SquareOfModulo() { return x * x + y * y; }

		template<typename T2> Point<T> operator += (Point<T2> & p) { x += p.x; y += p.y; return *this; }
		template<typename T2> Point<T> operator -= (Point<T2> & p) { x -= p.x; y -= p.y; return *this; }

		Point<T> operator += (Point<T> p) { x += p.x; y += p.y; return *this; }
		Point<T> operator -= (Point<T> p) { x -= p.x; y -= p.y; return *this; }

		Point<T> operator *= (T p) { x *= p; y *= p; return *this; }
		Point<T> operator /= (T p) { x /= p; y /= p; return *this; }

		Point<T> operator + () { return Point<T>(*this); }
		Point<T> operator - () { return Point<T>(-x, -y); }

		Point<T> Offset(T dx, T dy) { return Point<T>(x + dx, y + dy); }
		Point<T> & InnerOffset(T dx, T dy) { x += dx; y += dy; return *this; }
	};

	template<typename T> Point<T> operator + (Point<T> & p1, Point<T> & p2) { return Point<T>(p1.x + p2.x, p1.y + p2.y); }
	template<typename T> Point<T> operator - (Point<T> & p1, Point<T> & p2) { return Point<T>(p1.x - p2.x, p1.y - p2.y); }

	template<typename T> Point<T> operator * (Point<T> & p1, T p2) { return Point<T>(p1.x * p2, p1.y * p2); }
	template<typename T> Point<T> operator * (T p1, Point<T> & p2) { return Point<T>(p1 * p2.x, p1 * p2.y); }

	//两个点的乘积（点积） return p1.x * p2.y + p1.y * p2.x
	template<typename T> T operator * (Point<T> p1, Point<T> & p2) { return p1.x * p1.x + p2.y * p2.y; }

	//两个点的乘积（差积） return p1.x * p2.y - p1.y * p2.x
	template<typename T> T operator ^ (Point<T> p1, Point<T> & p2) { return p1.x * p2.y - p1.y * p2.x; }

	template<typename T> Point<T> operator / (Point<T> & p1, T p2) { return Point<T>(p1.x / p2, p1.y / p2); }

	template<typename T> bool operator == (Point<T> & p1, Point<T> & p2) { return p1.x == p2.x && p1.y == p2.y; }

	template<typename T> bool operator != (Point<T> & p1, Point<T> & p2) { return p1.x != p2.x || p1.y != p2.y; }

	template<typename T> ostream & operator << (ostream & o, Point<T> & p)
	{
		return o << p.x << " " << p.y;
	}

	template<typename T> istream & operator >> (istream & i, Point<T> & p)
	{
		return i >> p.x >> p.y;
	}

	template<typename T> T Abs(xmv::Point<T> & p) { return p.Modulo(); }

	template<typename T, typename T2> xmv::Point<T> Pow(xmv::Point<T> & x, T2 y)
	{
		T R = Pow(x.Modulo(), y);
		double alpha = std::atan2(x.y, x.x);
		return xmv::Point<T>(R * cos(alpha), R * sin(alpha));
	}

	template<typename T> xmv::Point<T> Sqrt(xmv::Point<T> & x) { return Pow(x, 0.5); }
}