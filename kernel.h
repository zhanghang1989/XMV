#pragma once

#include "matrix.h"

//创建方形平均滤波窗口
namespace xmv
{
	template<class T> Matrix<T> CreateAverageSquareKernel(int size)
	{
		return Matrix<T>(size, size, 1) / (size * size);
	}

	template<class T> Matrix<T> CreateAverageCircularKernel(int size)
	{
		Matrix<T> m = Matrix<T>(size, size);
		T * ptr;
		float cx = (size - 1) / 2.0f;
		float cy = cx;
		float r = cx + 0.5;
		float r2 = r * r;
		int count = 0;
		float x;
		float y;
		enumImage(m.Buffer(), size, size, 0, size, 0, size, ptr)
		{
			x = x_____ - cx;
			y = y_____ - cy;
			if (x * x + y * y < r2)
			{
				*ptr = 1;
				count++;
			}
			else
				*ptr = 0;
		}
		enumImageEnd;

		m /= count;

		return m;
	}

	template<class T> Matrix<T> CreateLaplacianKernel_141()
	{
		T k[9] = {0, 1, 0, 1, -4, 1, 0, 1, 0};
		return Matrix<T> (3, 3, k, 9);
	}

	template<class T> Matrix<T> CreateLaplacianKernel_181()
	{
		Matrix<T> k(3, 3, 1);
		k[1][1] = -8;
		return k;
	}

	template<class T> Matrix<T> CreateGaussianFilterKernel(int size, T delta = 0)
	{
		if (delta == 0) delta = Sqrt((size * size - 1) / 12.0);
		double k = -1.0 / (delta * delta * 2);
		//cout << k << endl;
		Matrix<T> m = Matrix<T>(size, size);
		T * ptr;
		float cx = (size - 1) / 2.0f;
		float cy = cx;
		float r = cx + 0.5;
		float r2 = r * r;

		float x, y;
		T sum = 0;
		T val;
		enumImage(m.Buffer(), size, size, 0, size, 0, size, ptr)
		{
			x = x_____ - cx;
			y = y_____ - cy;
			*ptr = val = exp((x * x + y * y) * k);
			sum += val;
		}
		enumImageEnd;


		m /= sum;

		return m;	}
	
}