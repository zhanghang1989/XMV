#pragma once
#include "point.h"

#define PI 3.141592653589793
namespace xmv
{
	template<class T> class Polygon
	{
	private:
		Point<T> * points;
		int capacity;
		int count;

	public:
		Polygon(int capacity = 1024)
		{
			this->capacity = capacity;
			this->count = 0;
			points = new Point<T>[capacity];
		}

		bool AddPoint(T x, T y)
		{
			if (count >= capacity) return false;
			points[count].x = x;
			points[count].y = y;
			++count;
		}

		bool AddPoint(Point<T> p)
		{
			if (count >= capacity) return false;
			points[count].x = p.x;
			points[count].y = p.y;
			++count;
		}

		Point<T> * Points() { return points; }

		void Clear(){ count = 0; }
		int Count() { return count; }

		Point<T> & operator[](int id)
		{
			return points[id];
		}

		bool In(T x, T y)
		{
			return In(Point<T>(x, y));
		}

		bool In(Point<T> pt)
		{
			double sa = 0;
			for (int i = 0; i < count; ++i)
			{
				Point<T> p1 = points[i];
				Point<T> p2 = points[(i + 1) % count];

				Point<T> c1 = p1 - pt;
				Point<T> c2 = p2 - pt;

				T dx = c1.x * c2.x + c1.y * c2.y;
				T dy = -c1.x * c2.y + c2.x * c1.y;

				double a = std::atan2(dy, dx);
				sa += a;
			}

			return fabs(sa) > PI;
		}

	};

	template <class T> istream & operator >> (istream & in, Polygon<T> & p)
	{
		p.Clear();
		int c = 0;
		in >> c;
		for (int i = 0; i < c; ++i)
		{
			Point<T> pt;
			in >> pt;
			p.AddPoint(pt);
		}
		return in;
	}

	template <class T> ostream & operator << (ostream & out, Polygon<T> & p)
	{
		out >> p.count;
		for (int i = 0; i < c; ++i)
			out >> p.points[i] >> endl;
		return out;
	}
}