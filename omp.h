#pragma once
#include "common.h"
#include "matrix.h"
#include "vector.h"

namespace xmv
{
	//function [A]=OMP(D,X,L); 
	//%=============================================
	//% Sparse coding of a group of signals based on a given 
	//% dictionary and specified number of atoms to use. 
	//% input arguments: 
	//%       D - the dictionary (its columns MUST be normalized).
	//%       X - the signals to represent
	//%       L - the Max. number of coefficients for each signal.
	//% output arguments: 
	//%       A - sparse coefficient matrix.
	//%=============================================
	template<typename T> void OMP(Matrix<T> & D, Matrix<T> & X, int L, Matrix<T> & A, T eps = 1E-10)
	{
		int n = X.Row();				//[n,P]=size(X);
		int P = X.Col();
		if (D.Row() != n) 
		{
			cout << "D.Row = " << D.Row() << endl;
			cout << "X.Row = " << X.Row() << endl;
			cout << "D.Row != n!!!\n";
			throw;
		}
		int K = D.Col();				//[n,K]=size(D);

		A = Matrix<T>(K, P, 0);

		int k = 0;
		for(k = 0; k < P; ++k)		//for k=1:1:P,
		{
			Vector<T> a;				//a=[];
			Vector<T> x = X.GetCol(k);	//x=X(:,k);

			if (x.IsZero())	continue;

			Vector<T> residual = x;		//residual=x;
			Vector<int> indx(L, 0);		//indx=zeros(L,1);
			int j = 0;
			for(j = 0; j < L; ++j)	//for j=1:1:L,
			{
				VectorH<T> proj = residual.ConjugateTranspose() * D;	// proj=D'*residual;
				T maxValue = 0;
				int maxValueIndex = 0;
				for(int ix = 0; ix < proj.Size(); ++ix)	// [maxVal,pos]=Max(Abs(proj));
				{
					if (Abs<T>(proj[ix]) > maxValue && ! indx.Contains(ix, 0, j)) 
					{
						maxValue = Abs<T>(proj[ix]);
						maxValueIndex = ix;				// pos=pos(1);
					}
				}
				indx[j] = maxValueIndex;				// indx(j)=pos;

				Matrix<T> DD = D.SelectColumns(indx.Cut(0, j + 1));
				////cout << DD << endl;
				//a = (!DD) * x;		//a=pinv(D(:,indx(1:j)))*x;

				a = xmv::SolveApproximateLinearEquations(DD, x);

				if (!a.InfinityCheck())
				{
					//cout << a << endl;
				}

				////cout << a << endl;
				residual = x - DD * a;						//residual=x-D(:,indx(1:j))*a;
				if (Sqrt(residual.SquareOfNorm2() / residual.Size()) < eps)				//if sum(residual.^2) < 1e-6
				{
					break;								//	break;
				}	//end
			}	//	end;
			Vector<T> temp(K, 0);		//temp=zeros(K,1);
			for(int ix = 0; ix < j; ++ix) temp[indx[ix]] = a[ix];		//temp(indx(1:j))=a;
			A.SetCol(k, temp);		//A(:,k)=sparse(temp);

			if (!temp.InfinityCheck())
			{
				//cout << temp << endl;
			}

			if (k % 100 == 0) 
				cout << "\r OMP ... " << k;


			//if (! (residual.Norm2() < 10))//cout << k << " " << residual.Norm2() << endl;
			//if (! (a.Norm2() < 100))//cout << k << " " << a << endl;
		}		//end;
		cout << "\n";
		//system("pause");
		//return;
	}


	//function [A]=OMP(D,X,L); 
	//%=============================================
	//% Sparse coding of a group of signals based on a given 
	//% dictionary and specified number of atoms to use. 
	//% input arguments: 
	//%       D - the dictionary (its columns MUST be normalized).
	//%       X - the signals to represent
	//%       L - the Max. number of coefficients for each signal.
	//% output arguments: 
	//%       return - sparse coefficient Vector.
	//%=============================================
	template<typename T> Vector<T> OMP(Matrix<T> & D, Vector<T> & X, int L, T eps = 1E-10)
	{
		int n = X.Size();				//[n,P]=size(X);
		if (D.Row() != n) throw;
		int K = D.Col();				//[n,K]=size(D);

		//A = Matrix<T>(K, P, 0);

		Vector<T> a;				//a=[];
		Vector<T> x = X;	//x=X(:,k);

		if (x.IsZero())	return Vector<T>(K, 0);

		Vector<T> residual = x;		//residual=x;
		Vector<int> indx(L, 0);		//indx=zeros(L,1);
		int j = 0;
		for(j = 0; j < L; ++j)	//for j=1:1:L,
		{
			VectorH<T> proj = residual.ConjugateTranspose() * D;	// proj=D'*residual;
			T maxValue = 0;
			int maxValueIndex = 0;
			for(int ix = 0; ix < proj.Size(); ++ix)	// [maxVal,pos]=Max(Abs(proj));
			{
				if (Abs<T>(proj[ix]) > maxValue) 
				{
					maxValue = proj[ix];
					maxValueIndex = ix;				// pos=pos(1);
				}
			}
			indx[j] = maxValueIndex;				// indx(j)=pos;

			Matrix<T> DD = D.SelectColumns(indx.Cut(0, j + 1));
			////cout << DD << endl;
			//a = (!DD) * x;		//a=pinv(D(:,indx(1:j)))*x;
			a = xmv::SolveApproximateLinearEquations(DD, x);

			////cout << a << endl;
			residual = x - DD * a;						//residual=x-D(:,indx(1:j))*a;
			if (residual.Norm2() / residual.Size() < eps)				//if sum(residual.^2) < 1e-6
			{
				break;								//	break;
			}	//end
		}	//	end;
		Vector<T> temp(K, 0);		//temp=zeros(K,1);
		for(int ix = 0; ix < j; ++ix) temp[indx[ix]] = a[ix];		//temp(indx(1:j))=a;

		return temp;
		//A.SetCol(k, temp);		//A(:,k)=sparse(temp);

		//system("pause");
		//return;
	}
}