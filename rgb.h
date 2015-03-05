#pragma once

#include <iostream>

namespace xmv
{
#pragma pack(push)
#pragma pack(3)

	template <typename T>
	struct BGR
	{
		T B;
		T G;
		T R;

		T Min() { return xmv::Min(xmv::Min(R, G), B); }
		T Max() { return xmv::Max(xmv::Max(R, G), B); }

		T Gray() { return (B + G + R) / 3; }

		BGR<T> & operator = (T val)	{ B = G = R = val; return *this; }
		template<typename T2> BGR<T> & operator = (BGR<T2> val)  { B = val.B; G = val.G; R = val.R; return *this; }
		BGR<T> & operator += (BGR<T> & v)	{ B += v.B; G += v.G; R += v.R; return *this; }
		BGR<T> & operator -= (BGR<T> & v)	{ B -= v.B; G -= v.G; R -= v.R; return *this; }
		BGR<T> & operator += (T & v)	{ B += v; G += v; R += v; return *this; }
		BGR<T> & operator -= (T & v)	{ B -= v; G -= v; R -= v; return *this; }
		BGR<T> & operator *= (T & v)	{ B *= v; G *= v; R *= v; return *this; }
		BGR<T> & operator /= (T & v)	{ B /= v; G /= v; R /= v; return *this; }
		template <typename T2> BGR<T> & operator += (T2 & v)	{ B += v; G += v; R += v; return *this; }
		template <typename T2> BGR<T> & operator -= (T2 & v)	{ B -= v; G -= v; R -= v; return *this; }
		template <typename T2> BGR<T> & operator *= (T2 & v)	{ B *= v; G *= v; R *= v; return *this; }
		template <typename T2> BGR<T> & operator /= (T2 & v)	{ B /= v; G /= v; R /= v; return *this; }

		BGR<T> operator ~ () { return BGR(~B, ~G, ~R); }

		BGR(T val) { B = G = R = val; }
		BGR(T b, T g, T r)
		{
			B = b; G = g; R = r;
		}
		BGR() {}
	};

#ifndef NORGB
	#undef RGB
	template <typename T> struct RGB
	{
		T R;
		T G;
		T B;

		T Min() { return xmv::Min(xmv::Min(R, G), B); }
		T Max() { return xmv::Max(xmv::Max(R, G), B); }

		RGB<T> & operator = (T val)  { B = G = R = val; return *this; }
		RGB<T> & operator += (RGB<T> & v)	{ B += v.B; G += v.G; R += v.R; return *this; }
		RGB<T> & operator -= (RGB<T> & v)	{ B -= v.B; G -= v.G; R -= v.R; return *this; }
		RGB<T> & operator += (T & v)	{ B += v; G += v; R += v; return *this; }
		RGB<T> & operator -= (T & v)	{ B -= v; G -= v; R -= v; return *this; }
		RGB<T> & operator *= (T & v)	{ B *= v; G *= v; R *= v; return *this; }
		RGB<T> & operator /= (T & v)	{ B /= v; G /= v; R /= v; return *this; }
		template <typename T2> RGB<T> & operator += (T2 & v)	{ B += v; G += v; R += v; return *this; }
		template <typename T2> RGB<T> & operator -= (T2 & v)	{ B -= v; G -= v; R -= v; return *this; }
		template <typename T2> RGB<T> & operator *= (T2 & v)	{ B *= v; G *= v; R *= v; return *this; }
		template <typename T2> RGB<T> & operator /= (T2 & v)	{ B /= v; G /= v; R /= v; return *this; }

		RGB<T> operator ~ () { return RGB(~R, ~G, ~B); }

		RGB(T val) { B = G = R = val; }
		RGB(T r, T g, T b) { B = b; G = g; R = r; }
		RGB() {}
	};

	template<typename T> inline RGB<T> Abs(RGB<T> x) 
	{
		RGB<T> v;
		v.R = Abs(x.R);
		v.G = Abs(x.G);
		v.B = Abs(x.B);
		return v;
	}
#endif

	template<typename T> inline BGR<T> Abs(BGR<T> x) 
	{
		BGR<T> v;
		v.R = Abs(x.R);
		v.G = Abs(x.G);
		v.B = Abs(x.B);
		return v;
	}

	template<typename T>
	ostream & operator << (ostream & o, BGR<T> & val)
	{
		return o << val.B << " " << val.G << " " << val.R;
	}

#ifndef NORGB
	template<typename T>
	ostream & operator << (ostream & o, RGB<T> & val)
	{
		return o << val.R << " " << val.G << " " << val.B;
	}
#endif

	template<typename T>
	istream & operator >> (istream & o, BGR<T> & val)
	{
		return o >> val.B >> val.G >> val.R;
	}

#ifndef NORGB
	template<typename T>
	istream & operator >> (istream & o, RGB<T> & val)
	{
		return o >> val.R >> val.G >> val.B;
	}


	template<typename T> bool operator <(RGB<T> v1, RGB<T> v2)
	{
		return (v1.R + v1.G + v1.B) < (v2.R + v2.G + v2.B);
	}

	template<typename T> bool operator >(RGB<T> v1, RGB<T> v2)
	{
		return (v1.R + v1.G + v1.B) > (v2.R + v2.G + v2.B);
	}

	template<typename T> bool operator >(RGB<T> v1, T v2)
	{
		return (v1.R + v1.G + v1.B) > v2 * 3;
	}

	template<typename T> bool operator < (RGB<T> v1, T v2)
	{
		return (v1.R + v1.G + v1.B) < v2 * 3;
	}
#endif

	template<typename T> bool operator >(BGR<T> v1, T v2)
	{
		return (v1.R + v1.G + v1.B) > v2 * 3;
	}

	template<typename T> bool operator < (BGR<T> v1, T v2)
	{
		return (v1.R + v1.G + v1.B) < v2 * 3;
	}

	template<typename T> bool operator <(BGR<T> v1, BGR<T> v2)
	{
		return (v1.R + v1.G + v1.B) < (v2.R + v2.G + v2.B);
	}

	template<typename T> bool operator >(BGR<T> v1, BGR<T> v2)
	{
		return (v1.R + v1.G + v1.B) > (v2.R + v2.G + v2.B);
	}

#ifndef NORGB
	template<typename T> bool operator <=(RGB<T> v1, RGB<T> v2)
	{
		return (v1.R + v1.G + v1.B) <= (v2.R + v2.G + v2.B);
	}

	template<typename T> bool operator >=(RGB<T> v1, RGB<T> v2)
	{
		return (v1.R + v1.G + v1.B) >= (v2.R + v2.G + v2.B);
	}
#endif

	template<typename T> bool operator <=(BGR<T> v1, BGR<T> v2)
	{
		return (v1.R + v1.G + v1.B) < (v2.R + v2.G + v2.B);
	}

	template<typename T> bool operator >=(BGR<T> v1, BGR<T> v2)
	{
		return (v1.R + v1.G + v1.B) >= (v2.R + v2.G + v2.B);
	}


	template<typename T> bool operator == (BGR<T> v1, T v2)
	{
		return v1.G == v2 && v1.R == v2 && v1.B == v2;
	}

#ifndef NORGB
	template<typename T> bool operator == (RGB<T> v1, T v2)
	{
		return v1.G == v2 && v1.R == v2 && v1.B == v2;
	}
#endif


	template<typename T> bool operator == (BGR<T> v1, BGR<T> v2)
	{
		return v1.G == v2.G && v1.R == v2.R && v1.B == v2.B;
	}

#ifndef NORGB
	template<typename T> bool operator == (RGB<T> v1, RGB<T> v2)
	{
		return v1.G == v2.G && v1.R == v2.R && v1.B == v2.B;
	}
#endif


	template<typename T> bool operator != (BGR<T> v1, T v2)
	{
		return !(v1 == v2);
	}

#ifndef NORGB
	template<typename T> bool operator != (RGB<T> v1, T v2)
	{
		return !(v1 == v2);
	}
#endif


	template<typename T> bool operator != (BGR<T> v1, BGR<T> v2)
	{
		return !(v1 == v2);
	}

#ifndef NORGB
	template<typename T> bool operator != (RGB<T> v1, RGB<T> v2)
	{
		return !(v1 == v2);
	}
#endif

	template<typename T> BGR<T> operator * (BGR<T> v1, T v2)
	{
		return BGR<T>(v1.B * v2, v1.G * v2, v1.R * v2);
	}

#ifndef NORGB
	template<typename T> RGB<T> operator * (RGB<T> v1, T v2)
	{
		return RGB<T>(v1.R * v2, v1.G * v2, v1.B * v2);
	}
#endif

	template<typename T> BGR<T> operator / (BGR<T> v1, T v2)
	{
		return BGR<T>(v1.B / v2, v1.G / v2, v1.R / v2);
	}

#ifndef NORGB
	template<typename T> RGB<T> operator / (RGB<T> v1, T v2)
	{
		return RGB<T>(v1.R / v2, v1.G / v2, v1.B / v2);
	}
#endif

#pragma pack(pop)

#pragma pack(push)
#pragma pack(4)
	template <typename T>
	struct ARGB
	{
		T A;
		T R;
		T G;
		T B;
	};
	template <typename T>
	struct BGRA
	{
		T B;
		T G;
		T R;
		T A;
	};
	template <typename T>
	struct ABGR
	{
		T A;
		T B;
		T G;
		T R;
	};
#pragma pack(pop)


}