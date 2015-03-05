#pragma once

#include "matrix.h"

//多项式拟合
namespace xmv
{
	//多项式拟合
	template<typename T> Vector<double> PolynomialFitting(T * xArray, T * yArray, int arraySize, int power)
	{
		Matrix<double> m(arraySize, power + 1);

		double * mPtr = m.Buffer() + (arraySize * (power + 1)) - 1;

		if (xArray)
		{
			T * xPtr = xArray + arraySize;
			double x = 1;
			int i, p;
			for(i = 0; i < arraySize; ++i, --xPtr)
			{
				x = 1;
				for(p = 0; p <= power; ++p, --mPtr, x *= *xPtr)
				{
					*mPtr = x;
				}
			}
		}
		else
		{
			double x = 1;
			int i, p;
			for(i = arraySize - 1; i >= 0; --i)
			{
				x = 1;
				for(p = 0; p <= power; ++p, --mPtr, x *= i)
				{
					*mPtr = x;
				}
			}
		}


		Vector<double> v(arraySize);
		double * vPtr = v.Buffer();
		double * vEnd = vPtr + arraySize;
		for(T * yPtr = yArray; vPtr < vEnd; ++vPtr, ++yPtr)
			*vPtr = *yPtr;

		return xmv::SolveApproximateLinearEquations(m, v);
	}

	//多项式拟合
	template<typename T> Vector<double> PolynomialFitting(T * yArray, int arraySize, int power)
	{
		return PolynomialFitting((T *)NULL, yArray, arraySize, power);
	}

	//二次函数拟合 y = ax^2 + bx + c
	template<typename T> void QuadraticFitting(T * xArray, T * yArray, int arraySize, double & a, double & b, double & c)
	{
		Vector<double> v = PolynomialFitting(xArray, yArray, arraySize, 2);
		double * vPtr = v.Buffer();
		a = *(vPtr++);
		b = *(vPtr++);
		c = *vPtr;
	}

	//直线拟合 y = ax + b
	template<typename T> void LinearFitting(T * yArray, int arraySize, double & a, double & b)
	{
		Vector<double> v = PolynomialFitting(yArray, arraySize, 1);
		a = v.GetValue(0);
		b = v.GetValue(1);
	}

	//直线拟合 y = ax + b
	template<typename T> void LinearFitting(T * xArray, T * yArray, int arraySize, double & a, double & b)
	{
		Vector<double> v = PolynomialFitting(xArray, yArray, arraySize, 1);
		a = v.GetValue(0);
		b = v.GetValue(1);
	}


	//二次函数拟合 y = ax^2 + bx + c
	template<typename T> void QuadraticFitting(T * yArray, int arraySize, double & a, double & b, double & c)
	{
		Vector<double> v = PolynomialFitting(yArray, arraySize, 2);
		double * vPtr = v.Buffer();
		a = *(vPtr++);
		b = *(vPtr++);
		c = *vPtr;
	}



}