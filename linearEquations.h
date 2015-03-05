#pragma once

#include "matrix.h"
#include "vector.h"

namespace xmv
{
	//解满秩线性方程
	template <typename T> Vector<T> SolveLinearEquations(Matrix<T> & m, Vector<T> & v)
	{
		Matrix<T> am = xmv::CreateAugmentedMatrix(m, v);

		am.GaussianElimination(true);

		Vector<T> rv = am.GetCol(am.Col() - 1);
		return rv;
	}

	//用最小二乘法获得过完备线性方程的最佳近似解
	template <typename T> Vector<T> SolveApproximateLinearEquations(Matrix<T> & m, Vector<T> & v)
	{
		Matrix<T> mT = m.Transpose();
		return SolveLinearEquations(mT * m, mT * v);
	}



}