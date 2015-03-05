#pragma once

#include <istream>
#include "math.h"
#include "convert.h"

namespace xmv
{
	template <typename T>
	class Complex
	{
	public: 
		T real;
		T image;

		Complex() { }
		Complex(T r, T i) { real = r; image = i; }
		Complex(T r) { real = r; image = 0;}

		T Modulo() { return Sqrt<T>(real * real + image * image); }
		T SquareOfModulo() { return real * real + image * image; }

		template<typename T2> Complex<T> operator += (Complex<T2> & c) { real += c.real; image += c.image; return *this; }
		template<typename T2> Complex<T> operator -= (Complex<T2> & c) { real -= c.real; image -= c.image; return *this; }
		template<typename T2> Complex<T> operator *= (Complex<T2> & c) { T tr = real * c.real - image * c.image; image = real * c.image + image * c.real; real = tr; return *this; }
		template<typename T2> Complex<T> operator /= (Complex<T2> & c) 
		{
			T r = c.SquareOfModulo();
			T tr = (real * c.real + image * c.image) / r;
			image = (- real * c.image + image * c.real) / r;
			real = tr; 
			return *this; 
		}
		Complex<T> operator += (T c) { real += c; return *this; }
		Complex<T> operator -= (T c) { real -= c; return *this; }
		Complex<T> operator *= (T c) { real *= c; real *= c; return *this; }
		Complex<T> operator /= (T c) { real /= c; real /= c; return *this; }
		Complex<T> operator + () { return Complex<T>(*this); }
		Complex<T> operator - () { return Complex<T>(-real, -image); }
		Complex<T> operator ~ () { return Complex<T>(real, -image); }
	};

	template<typename T> Complex<T> operator + (Complex<T> & c1, Complex<T> & c2) { return Complex<T>(c1.real + c2.real, c1.image + c2.image); }
	template<typename T> Complex<T> operator + (Complex<T> & c1, T c2) { return Complex<T>(c1.real + c2, c1.image); }
	template<typename T> Complex<T> operator + (T c1, Complex<T> & c2) { return Complex<T>(c1 + c2.real, c2.image); }

	template<typename T> Complex<T> operator - (Complex<T> & c1, Complex<T> & c2) { return Complex<T>(c1.real - c2.real, c1.image - c2.image); }
	template<typename T> Complex<T> operator - (Complex<T> & c1, T c2) { return Complex<T>(c1.real - c2, c1.image); }
	template<typename T> Complex<T> operator - (T c1, Complex<T> & c2) { return Complex<T>(c1 - c2.real, -c2.image); }

	template<typename T> Complex<T> operator * (Complex<T> & c1, Complex<T> & c2) { return Complex<T>(c1.real * c2.real - c1.image * c2.image, c1.real * c2.image + c2.real * c1.image); }
	template<typename T> Complex<T> operator * (Complex<T> & c1, T c2) { return Complex<T>(c1.real * c2, c1.image * c2); }
	template<typename T> Complex<T> operator * (T c1, Complex<T> & c2) { return Complex<T>(c1 * c2.real, c1 * c2.image); }

	template<typename T> Complex<T> operator / (Complex<T> & c1, Complex<T> & c2) 
	{
		T r = c2.SquareOfModulo();
		return Complex<T>((c1.real * c2.real + c1.image * c2.image) / r, (- c1.real * c2.image + c2.real * c1.image) / r); 
	}
	template<typename T> Complex<T> operator / (Complex<T> & c1, T c2) { return Complex<T>(c1.real / c2, c1.image / c2); }
	template<typename T> Complex<T> operator / (T c1, Complex<T> & c2) 
	{ 
		T r = c2.SquareOfModulo();
		return Complex<T>(c1 * c2.real / r, - c1.real * c2.image / r);
	}

	template<typename T> bool operator == (Complex<T> & c1, Complex<T> & c2) { return c1.real ==c2.real && c1.image == c2.image; }
	template<typename T> bool operator == (Complex<T> & c1, T c2) { return c1.real == c2 && c1.image == 0; }
	template<typename T> bool operator == (T c1, Complex<T> & c2) { return c1 == c2.real && c2.image == 0; }

	template<typename T> bool operator != (Complex<T> & c1, Complex<T> & c2) { return c1.real != c2.real || c1.image != c2.image; }
	template<typename T> bool operator != (Complex<T> & c1, T c2) { return c1.real != c2 || c1.image != 0; }
	template<typename T> bool operator != (T c1, Complex<T> & c2) { return c1 != c2.real || c2.image != 0; }

	template<typename T> ostream & operator << (ostream & o, Complex<T> & c) 
	{ 
		o << c.real;
		if (c.image == 0) return o;
		else if (c.image > 0) o << "+";
		o << c.image << "i"; 
		return o;
	}

	template<typename T> Complex<T> Conjugate(Complex<T> & x) { return ~x; }
	template<typename T> T Conjugate(T & x) { return x; }

	template<typename T> int Sign(Complex<T> v) { if (v == 0) return 0; else return 1; }

	template<typename T> inline T Abs(xmv::Complex<T> & c) { return c.Modulo(); }

	template<typename T, typename T2> xmv::Complex<T> Pow(xmv::Complex<T> & x, T2 y)
	{
		T R = std::Pow(x.Modulo(), y);
		double alpha = std::atan2(x.image, x.real);
		return xmv::Complex<T>(R * cos(alpha), R * sin(alpha));
	}

	template<typename T> xmv::Complex<T> Sqrt(xmv::Complex<T> & x) { return Pow<Complex<T>>(x, 0.5); }
}