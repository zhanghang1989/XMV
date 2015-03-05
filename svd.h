/*
*
* This program is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 2 or any later version.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
* more details. A copy of the GNU General Public License is available at:
* http://www.fsf.org/licensing/licenses
*/


/*****************************************************************************
*                                    svd.h
*
* Class template of Singular Value Decomposition.
*
* For an m-by-n matrix A, the singular value decomposition is an m-by-m
* orthogonal matrix U, an m-by-n diagonal matrix S, and an n-by-n orthogonal
* matrix V so that A = U*S*V^T.
*
* The singular values, sigma[k] = S[k][k], are ordered so that sigma[0] >=
* sigma[1] >= ... >= sigma[n-1].
*
* For economy size, denotes p = Min(m,n), then U is m-by-p, S is p-by-p,
* and V is n-by-p, this file provides the economic decomposition format.
*
* The singular value decompostion always exists, so the constructor will
* never fail. The matrix condition number and the effective numerical rank
* can be computed from this decomposition.
*
* Adapted form Template Numerical Toolkit.
*
*****************************************************************************/


#pragma once
#include <cmath>
#include "matrix.h"
#include "common.h"

namespace xmv
{
	template <typename T> class Matrix;
	template <typename T> class Vector;
	template <typename T>
	class SVD
	{
	public:

		SVD();
		~SVD();

		void dec(Matrix<T> &A);
		Matrix<T> getU();
		Matrix<T> getV();
		Matrix<T> getSM();
		Vector<T> getSV();

		T norm2();
		T cond();
		int  rank();
		T eps;

	private:

		// the orthogonal matrix and singular value vector
		Matrix<T> U;
		Matrix<T> V;
		Vector<T> S;

		// docomposition for matrix with rows >= columns
		void decomposition(Matrix<T>&, Matrix<T>&,
			Vector<T>&, Matrix<T>&);

	};
	// class SVD

	/**
	* constructor and destructor
	*/
	template<typename T>
	SVD<T>::SVD()
	{
		eps = 1E-10;
	}

	template<typename T>
	SVD<T>::~SVD()
	{
	}


	/**
	* Making singular decomposition.
	*/
	template <typename T>
	void SVD<T>::dec(Matrix<T> & A)
	{
		int m = A.Row(),
			n = A.Col(),
			p = Min(m, n);

		U = Matrix<T>(m, p, 0);
		V = Matrix<T>(n, p, 0);
		S = Vector<T>(p, 0);
		if (m >= n)
		{
			Matrix<T> B(A);
			decomposition(B, U, S, V);
		}
		else
		{
			Matrix<T> B(A.ConjugateTranspose());
			decomposition(B, V, S, U);
		}
	}


	/**
	* Making singular decomposition of m >= n.
	*/
	template <typename T>
	void SVD<T>::decomposition(Matrix<T> &B, Matrix<T> &U, Vector<T> &S, Matrix<T> &V)
	{
		int m = B.Row(),
			n = B.Col();

		T * b = B.Buffer();
		T * u = U.Buffer();
		T * s = S.Buffer();
		T * v = V.Buffer();

		T t;

		Vector<T> e(n);
		Vector<T> work(m);

		// boolean
		int wantu = 1;
		int wantv = 1;

		// Reduce A to bidiagonal form, storing the diagonal elements
		// in s and the super-diagonal elements in e.
		int nct = Min(m - 1, n);
		int nrt = Max(0, n - 2);
		int i = 0,
			j = 0,
			k = 0;

		int maxNCRT = Max(nct, nrt);
		T * pbkk = b;
		T * pvkk = v;
		T * pukk = u;
		T * pbnn;
		T * pvij;
		T * pvik;
		T * pvkj;
		T * pbij;
		T * pbik;
		T * pbkj;
		T * puij;
		T * pvij2;
		T * puik;
		T * pukj;
		T * pvk1k;
		T * puij2;
		T * pe;
		T * pe0 = e.Buffer();
		T * pek = pe0;
		T * peka1 = pe0 + 1;
		T * pei;
		T * pej;
		T * pw0 = work.Buffer();
		T * pw;
		//T & snct = s[nct];
		//T & sn_1 = s[n - 1];
		//T & bnctnct = b[nct * n + nct];
		//T & enrt = pe0[nrt];
		//T & bnrtn_1 = b[nrt * n + n - 1];
		//T & en_1 = pe0[n - 1];
		T * snct = s + nct;
		T * sn_1 = s + n - 1;
		T * bnctnct = b + nct * n + nct;
		T * enrt = pe0 + nrt;
		T * bnrtn_1 = b + nrt * n + n - 1;
		T * en_1 = pe0 + n - 1;
		T * psk = s; // s_k_
		T * psk2;
		T scale;
	//	T t;

		T sp;
		T spm1;
		T epm1;
		T sk;
		T ek;
		T b_;
		T c;
		T shift;
		T f;
		T g;

		for (k = 0; k < maxNCRT; ++k, pbkk += n + 1, ++pek, ++peka1, pvkk += n + 1, pukk += n + 1, ++psk)
		{
			T & s_k_ = *psk;
			//T &	bkk		= *pbkk;
			T & ek = *pek;

			if (k < nct)
			{
				// Compute the transformation for the k-th column and
				// place the k-th diagonal in s[k].
				// Compute 2-norm of k-th column without under/overflow.
				s_k_ = 0;
				pbik = pbkk;
				for (i = k; i < m; ++i, pbik += n)
				{
					s_k_ = hypot(s_k_, *pbik /*B[i][k]*/);
				}

				if (s_k_ != 0)
				{
					if (*pbkk/*B[k][k]*/ < 0)
						s_k_ = -s_k_;

					pbik = pbkk;
					for (i = k; i < m; ++i, pbik += n)
						*pbik /= s_k_;
					++(*pbkk)/*B[k][k] += 1*/;
				}
				s_k_ = -s_k_;
			}

			pbkj = pbkk + 1;
			pe = pe0 + k + 1;
			for (j = k + 1; j < n; ++j, ++pbkj, ++pe)
			{
				if ((k < nct) && (s_k_ != 0))
				{
					// apply the transformation
					t = 0;
					pbik = pbkk;
					pbij = pbkk - k + j;
					for (i = k; i < m; ++i, pbik += n, pbij += n)
						t += *pbik * *pbij; // B[i][j];

					pbik = pbkk;
					t = -t / *pbkk;
					pbij = pbkk - k + j;
					for (i = k; i < m; ++i, pbik += n, pbij += n)
						*pbij /*B[i][j]*/ += t * *pbik;
				}
				*pe /*e[j]*/ = *pbkj;//  B[k][j];
			}

			// Place the transformation in U for subsequent back
			// multiplication.
			if (wantu & (k < nct))
			{
				pbik = pbkk;
				puik = pukk;
				for (i = k; i < m; ++i, pbik += n, puik += n)
					*puik /*U[i][k]*/ = *pbik;// B[i][k];
			}

			if (k < nrt)
			{
				// Compute the k-th row transformation and place the
				// k-th super-diagonal in e[k].
				// Compute 2-norm without under/overflow.
				ek = 0;
				pei = pe0 + k + 1;
				//T * ei2 = ei;
				for (i = k + 1; i < n; ++i, ++pei)
					ek = hypot(ek, *pei);

				if (ek)
				{
					if (/*e[k + 1]*/*peka1 < 0)
						ek = -ek;

					pei = pe0 + k + 1;
					for (i = k + 1; i < n; ++i, ++pei)
						*pei /= ek;
					++*peka1;// += 1;
				}
				ek = -ek;

				if ((k + 1 < m) && (ek))
				{
					// apply the transformation
					pw = pw0 + k + 1;
					for (i = k + 1; i < m; ++i, ++pw)
						*pw/* work[i]*/ = 0;

					pej = pe0 + k + 1;
					pbkj = pbkk + n + 1;
					for (j = k + 1; j < n; ++j, ++pej, ++pbkj)
					{
						pbij = pbkj;
						pw = pw0 + k + 1;
						for (i = k + 1; i < m; ++i, pbij += n, ++pw)
							*pw += *pej * *pbij;// B[i][j];
					}

					pej = pe0 + k + 1;
					pbkj = pbkk + n + 1;
					for (j = k + 1; j < n; ++j, ++pej, ++pbkj)
					{
						pbij = pbkj;
						t = -*pej / *peka1;
						pw = pw0 + k + 1;
						for (i = k + 1; i < m; ++i, pbij += n, ++pw)
							*pbij += t * *pw;
					}
				}

				// Place the transformation in V for subsequent
				// back multiplication.
				pei = pe0 + k + 1;
				if (wantv)
				{
					pvik = pvkk + n;
					for (i = k + 1; i < n; ++i, ++pei, pvik += n)
						*pvik = *pei;
				}
			}
		}

		// Set up the final bidiagonal matrix or order p.
		int p = n;

		if (nct < n)
			*snct = *bnctnct;
		if (m < p)
			*sn_1 = 0;
		/*s[p - 1] = 0;*/

		if (nrt + 1 < p)
			*enrt = *bnrtn_1;
		//e[nrt] = B[nrt][p - 1];
		//e[p - 1] = 0;
		*en_1 = 0;

		// if required, generate U
		if (wantu)
		{
			puik = u + nct;
			for (i = 0; i < m; ++i, puik += n)
			{
				puij = puik;
				for (j = nct; j < n; ++j, ++puij)
					*puij = i == j ? 1 : 0;
				//U[i][j] = 0;
				//U[j][j] = 1;
				//*(puik + j) = 1;
			}

			psk = s + nct - 1;
			pukk = u + (nct - 1) * n + nct - 1;
			for (k = nct - 1; k >= 0; --k, --psk, pukk -= n + 1)
			{
				if (*psk)
				{
					for (j = k + 1; j < n; ++j)
					{
						t = 0;
						puij = pukk + j - k;
						puik = pukk;
						for (i = k; i < m; ++i, puik += n, puij += n)
							t += *puik /*U[i][k]*/ * *puij;// U[i][j];
						t = -t / *pukk;//U[k][k];

						puik = pukk;
						puij = pukk + j - k;
						for (i = k; i < m; ++i, puik += n, puij += n)
							*puij/*U[i][j]*/ += t * *puik;// U[i][k];
					}

					puik = pukk;
					for (i = k; i < m; ++i, puik += n)
						*puik = -*puik;//U[i][k] = -U[i][k];

					++*pukk;//U[k][k] = 1 + U[k][k];

					puik = u + k;
					for (i = 0; i < k - 1; ++i, puik += n)
						*puik = 0; // U[i][k] = 0;
				}
				else
				{
					puik = u + k;
					for (i = 0; i < m; ++i, puik += n)
						*puik = 0;
					//U[i][k] = 0;
					*pukk = 1;
					//U[k][k] = 1;
				}
			}
		}

		// if required, generate V
		if (wantv)
		{
			pe = pe0 + n - 1;
			pvkk = v + (n - 1) * n + n - 1;
			pvk1k = pvkk + n;
			for (k = n - 1; k >= 0; --k, --pe, pvkk -= n + 1, pvk1k -= n + 1)
			{
				if ((k < nrt) && (*pe/*e[k]*/))
				{
					pvkj = pvkk + 1;
					for (j = k + 1; j < n; ++j, ++pvkj)
					{
						t = 0;
						pvik = pvkk + n;
						pvij = pvkj + n;
						for (i = k + 1; i < n; ++i, pvik += n, pvij += n)
							t += *pvik /*V[i][k] */* *pvij /*V[i][j]*/;
						t = -t / *pvk1k;//V[k + 1][k];

						pvik = pvkk + n;
						pvij = pvkj + n;
						for (i = k + 1; i < n; ++i, pvik += n, pvij += n)
							*pvij/*V[i][j]*/ += t * *pvik/*V[i][k]*/;
					}
				}

				pvik = v + k;
				for (i = 0; i < n; ++i, pvik += n)
					*pvik = 0;
				/*V[i][k] = 0;*/
				//V[k][k] = 1;
				*pvkk = 1;
			}
		}

		// main iteration loop for the singular values
		int pp = p - 1;
		int iter = 0;
		//double eps = Pow<T>( 2.0, -52.0 );

		while (p > 0)
		{
			int k = 0;
			int kase = 0;

			// Here is where a test for too many iterations would go.
			// This section of the program inspects for negligible
			// elements in the s and e arrays. On completion the
			// variables kase and k are set as follows.
			// kase = 1     if s(p) and e[k-1] are negligible and k<p
			// kase = 2     if s(k) is negligible and k<p
			// kase = 3     if e[k-1] is negligible, k<p, and
			//				s(k), ..., s(p) are not negligible
			// kase = 4     if e(p-1) is negligible (convergence).

			pek = pe0 + p - 2;
			psk = s + p - 2;
			psk2 = psk + 1;
			for (k = p - 2; k >= -1; --k, --pek, --psk, --psk2)
			{
				if (k == -1)
					break;

				if (Abs<T>(*pek) <= eps*(Abs<T>(*psk/*s[k]*/) + Abs<T>(*psk2/*s[k + 1]*/)))
				{
					*pek = 0;
					break;
				}
			}

			if (k == p - 2)
				kase = 4;
			else
			{
				int ks;
				pek = pe0 + p - 1;
				psk = s + p - 1;
				psk2 = psk - 1;
				for (ks = p - 1; ks >= k; --ks, --pek, --psk, --psk2)
				{
					if (ks == k)
						break;

					t = ((ks != p) ? Abs<T>(*pek) : 0) +
						((ks != k + 1) ? Abs<T>(*psk2) : 0);

					if (Abs<T>(*psk) <= eps*t)
					{
						s[ks] = 0;
						break;
					}
				}

				if (ks == k)
					kase = 3;
				else if (ks == p - 1)
					kase = 1;
				else
				{
					kase = 2;
					k = ks;
				}
			}
			k++;

			// Perform the task indicated by kase.
			switch (kase)
			{
				// deflate negligible s(p)
			case 1:
			{
					  T f = e[p - 2];
					  e[p - 2] = 0;

					  for (j = p - 2; j >= k; --j)
					  {
						  t = hypot(s[j], f);
						  T cs = s[j] / t;
						  T sn = f / t;
						  s[j] = t;

						  if (j != k)
						  {
							  f = -sn * e[j - 1];
							  e[j - 1] = cs * e[j - 1];
						  }

						  if (wantv)
						  for (i = 0; i < n; ++i)
						  {
							  t = cs*V[i][j] + sn*V[i][p - 1];
							  V[i][p - 1] = -sn*V[i][j] + cs*V[i][p - 1];
							  V[i][j] = t;
						  }
					  }
			}
				break;

				// split at negligible s(k)
			case 2:
			{
					  T f = e[k - 1];
					  e[k - 1] = 0;

					  for (j = k; j < p; ++j)
					  {
						  t = hypot(s[j], f);
						  T cs = s[j] / t;
						  T sn = f / t;
						  s[j] = t;
						  f = -sn * e[j];
						  e[j] = cs * e[j];

						  if (wantu)
						  for (i = 0; i < m; ++i)
						  {
							  t = cs*U[i][j] + sn*U[i][k - 1];
							  U[i][k - 1] = -sn*U[i][j] + cs*U[i][k - 1];
							  U[i][j] = t;
						  }
					  }
			}
				break;

				// perform one qr step
			case 3:
			{
					  // calculate the shift
					  scale = Max(Max(Max(Max(
						  Abs<T>(s[p - 1]), Abs<T>(s[p - 2])), Abs<T>(e[p - 2])),
						  Abs<T>(s[k])), Abs<T>(e[k]));
					  sp   =    s[p - 1] / scale;
					  spm1 = s[p - 2] / scale;
					  epm1 = e[p - 2] / scale;
					  sk   =   s[k] / scale;
					  ek   =   e[k] / scale;
					  b_    =   ((spm1 + sp)*(spm1 - sp) + epm1*epm1) / 2.0;
					  c    =   (sp*epm1) * (sp*epm1);
					  shift = 0;

					  if ((b_ != 0) || (c != 0))
					  {
						  shift = Sqrt<T>(b_*b_ + c);
						  if (b_ < 0)
							  shift = -shift;
						  shift = c / (b_ + shift);
					  }
					  f = (sk + sp)*(sk - sp) + shift;
					  g = sk * ek;

					  // chase zeros
					  pvik = v + k;
					  puik = u + k;
					  for (j = k; j < p - 1; ++j, pvik++, puik++)
					  {
						   t = hypot(f, g);
						  T cs = f / t;
						  T sn = g / t;
						  if (j != k)
							  e[j - 1] = t;

						  f = cs*s[j] + sn*pe0[j];
						  e[j] = cs*e[j] - sn*s[j];
						  g = sn * s[j + 1];
						  s[j + 1] = cs * s[j + 1];

						  if (wantv)
						  {
							  pvij = pvik;
							  pvij2 = pvij + 1;
							  for (i = 0; i < n; ++i, pvij += n, pvij2 += n)
							  {
								  t = cs* *pvij + sn * *pvij2/*V[i][j + 1]*/;
								  *pvij2 /*V[i][j + 1]*/ = -sn * *pvij + cs * *pvij2/*V[i][j + 1]*/;
								  *pvij = t;
							  }
						  }

						  t = hypot(f, g);
						  cs = f / t;
						  sn = g / t;
						  s[j] = t;
						  f = cs*e[j] + sn*s[j + 1];
						  s[j + 1] = -sn*e[j] + cs*s[j + 1];
						  g = sn * e[j + 1];
						  e[j + 1] = cs * e[j + 1];

						  if (wantu && (j < m - 1))
						  {
							  puij = puik;
							  puij2 = puij + 1;
							  for (i = 0; i < m; ++i, puij += n, puij2 += n)
							  {
								  t = cs * *puij + sn * *puij2;
								  //t = cs*U[i][j] + sn*U[i][j + 1];
								  *puij2 = -sn * *puij + cs * *puij2;
								  // U[i][j + 1] = -sn*U[i][j] + cs*U[i][j + 1];
								  *puij = t;
								  //  U[i][j] = t;*/
							  }
						  }
					  }
					  e[p - 2] = f;
					  iter = iter + 1;
			}
				break;

				// convergence
			case 4:
			{
					  // Make the singular values positive.
					  if (s[k] <= 0)
					  {
						  s[k] = (s[k] < 0) ? -s[k] : 0;
						  if (wantv)
						  {
							  pvik = v + k;
							  for (i = 0; i <= pp; ++i, pvik += n)
								  *pvik = -*pvik;
						  }
					  }

					  // Order the singular values.
					  while (k < pp)
					  {
						  if (s[k] >= s[k + 1])
							  break;

						  t = s[k];
						  s[k] = s[k + 1];
						  s[k + 1] = t;

						  if (wantv && (k < n - 1))
						  {
							  pvik = v + k;
							  pvij = v + k + 1;
							  for (i = 0; i < n; ++i, pvik += n, pvij += n)
								  swap(*pvik, *pvij); // V[i][k], V[i][k + 1]);
						  }

						  if (wantu && (k < m - 1))
						  {
							  pvik = v + k;
							  pvij = v + k + 1;
							  for (i = 0; i < m; ++i)
								  swap(*pvik, *pvij);
							  //swap(U[i][k], U[i][k + 1]);
						  }
						  k++;
					  }
					  iter = 0;
					  p--;
			}
				break;
			}
		}
	}


	/**
	* Get the left singular vectors.
	*/
	template<typename T>
	inline Matrix<T> SVD<T>::getU()
	{
		return U;
	}


	/**
	* Get the singular values matrix.
	*/
	template<typename T>
	inline Matrix<T> SVD<T>::getSM()
	{
		int N = S.Size();
		Matrix<T> tmp(N, N, 0);
		for (int i = 0; i < N; ++i)
			tmp[i][i] = S[i];

		return tmp;
	}


	/**
	* Get the singular values vector.
	*/
	template<typename T>
	inline Vector<T> SVD<T>::getSV()
	{
		return S;
	}


	/**
	* Get the right singular vectors.
	*/
	template<typename T>
	inline Matrix<T> SVD<T>::getV()
	{
		return V;
	}


	/**
	* Two norm (Max(S)).
	*/
	template <typename T>
	inline T SVD<T>::norm2()
	{
		return S[0];
	}


	/**
	* Two norm of condition number (Max(S)/Min(S)).
	*/
	template <typename T>
	inline T SVD<T>::cond()
	{
		return (S[0] / S(S.Size()));
	}


	/**
	* Effective numerical matrix rank.
	*/
	template <typename T>
	int SVD<T>::rank()
	{
		int N = S.Size();
		double tol = N * S[0] * eps;
		int r = 0;

		for (int i = 0; i<N; ++i)
		if (S[i] > tol)
			r++;

		return r;
	}



	/*****************************************************************************
	*                               svd_test.cpp
	*
	* SVD class testing.
	*
	* Zhang Ming, 2010-01 (revised 2010-08), Xi'an Jiaotong University.
	*****************************************************************************/

	//
	//int main()
	//{
	//	//    Matrix<Type> A(4,4);
	//	//	A(1,1) = 16;    A(1,2) = 2;     A(1,3) = 3;     A(1,4) = 13;
	//	//	A(2,1) = 5;     A(2,2) = 11;    A(2,3) = 10;    A(2,4) = 8;
	//	//	A(3,1) = 9;     A(3,2) = 7;     A(3,3) = 6;     A(3,4) = 12;
	//	//	A(4,1) = 4;     A(4,2) = 14;    A(4,3) = 15;    A(4,4) = 1;
	//
	//	//    Matrix<Type> A(4,2);
	//	//	A(1,1) = 1;     A(1,2) = 2;
	//	//	A(2,1) = 3;     A(2,2) = 4;
	//	//	A(3,1) = 5;     A(3,2) = 6;
	//	//	A(4,1) = 7;     A(4,2) = 8;
	//
	//	Matrix<Type> A(2,4);
	//	A(1,1) = 1;     A(1,2) = 3;     A(1,3) = 5;     A(1,4) = 7;
	//	A(2,1) = 2;     A(2,2) = 4;     A(2,3) = 6;     A(2,4) = 8;
	//
	//	SVD<Type> svd;
	//	svd.dec(A);
	//
	//	Matrix<Type> U = svd.getU();
	//	Matrix<Type> V = svd.getV();
	//	Matrix<Type> S = svd.getSM();
	//
	//	//cout << setiosflags(ios::fixed) << setprecision(4);
	//	//cout << "Matrix--A: " << A << endl;
	//	//cout << "Matrix--U: " << U << endl;
	//	//cout << "Vector--S: " << S << endl;
	//	//cout << "Matrix--V: " << V << endl;
	//	//cout << "Matrix--A - U * S * V^T:  "
	//		<< A- U*S*trT(V) << endl;
	//	//         << A- U*multTr(S,V) << endl;
	//
	//	//cout << "The rank of A : " << svd.rank() << endl << endl;
	//	//cout << "The condition number of A : " << svd.cond() << endl << endl;
	//
	//	return 0;
	//}
	//运行结果：
	//
	//
	//	Matrix--A: size: 2 by 4
	//	1.0000  3.0000  5.0000  7.0000
	//	2.0000  4.0000  6.0000  8.0000
	//
	//	Matrix--U: size: 2 by 2
	//	0.6414  -0.7672
	//	0.7672  0.6414
	//
	//	Vector--S: size: 2 by 2
	//	14.2691 0.0000
	//	0.0000  0.6268
	//
	//	Matrix--V: size: 4 by 2
	//	0.1525  0.8226
	//	0.3499  0.4214
	//	0.5474  0.0201
	//	0.7448  -0.3812
	//
	//	Matrix--A - U * S * V^T:  size: 2 by 4
	//	-0.0000 0.0000  -0.0000 0.0000
	//	-0.0000 -0.0000 -0.0000 -0.0000
	//
	//	The rank of A : 2
	//
	//	The condition number of A : 22.7640
	//
	//
	//	Process returned 0 (0x0)   execution time : 0.094 s
	//	Press any key to continue.
	//
}