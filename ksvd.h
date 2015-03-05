#include <list>
#include <cstdlib>
#include "vector.h"
#include "matrix.h"
#pragma once

#include "common.h"
#include "memory.h"



namespace xmv
{
	typedef double real;
	///KSVD算法类
	template <typename T>
	class KSVD
	{
	private:
		list<Vector<T> *> samples;		//样本

		Matrix<T> * dictionaryPtr;		//字典矩阵：n x K，n为样本维度，K为码本数量（即最后要确定的码本）
		Matrix<T> * samplesPtr;			//样本矩阵：n x N, n为样本维度，N为样本数量
		Matrix<T> * coefficientPtr;		//系数矩阵：K x N, K为码本数量，N为样本数量
		int codeBookCount;				//K，码本数量
		int sampleCount;				//N，样本数量
		int sampleDimension;			//n, 样本维度

		ClusterAnalyzer<T> ca;
	public:
		Matrix<T> Dictionary()
		{
			return *dictionaryPtr;
		}


	public:
		KSVD()				//构造函数
			: dictionaryPtr(NULL), samplesPtr(NULL), coefficientPtr(NULL),
			codeBookCount(0), sampleCount(0), sampleDimension(0)
		{
		}

		//析构函数
		~KSVD()
		{
			//DELETE samples
			ClearSamples();

			if (dictionaryPtr)	DELETE dictionaryPtr;
			if (samplesPtr)		DELETE samplesPtr;
			if (coefficientPtr)	DELETE coefficientPtr;
		}

	public:
		//返回码本数量
		int CodebookCount() { return this->codeBookCount; }

		//通过一个向量添加一个样本
		void AddSample(Vector<T> & sample)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(Vector<T> & sample)") Vector<T>(sample));
			ca.AddSample(sample);
		}

		//通过一个数组来添加一个样本。
		void AddSample(T * sampleArray, int arraySize)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(T * sampleArray, int arraySize)")
				Vector<T>(sampleArray, arraySize));
			ca.AddSample(sampleArray, arraySize);
		}

		//通过一个向量添加一个样本
		template<typename ST>
		void AddSample(Vector<ST> & sample)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(Vector<ST> & sample)") Vector<T>(sample));
			ca.AddSample(sample);
		}

		//通过一个数组来添加一个样本。
		template<typename ST>
		void AddSample(ST * sampleArray, int arraySize)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(ST * sampleArray, int arraySize)")
				Vector<T>(sampleArray, arraySize));
			ca.AddSample(sampleArray, arraySize);
		}

		//删除所有的样本
		void ClearSamples()
		{
			list<Vector<T> *>::iterator i = samples.begin();
			list<Vector<T> *>::iterator end = samples.end();

			for (; i != end; ++i)
				DELETE(*i);

			samples.clear();
			ca.ClearSamples();
		}

	private:
		//创建样本数组
		bool CreateSamplesMatrix()
		{
			//删除原有样本
			if (samplesPtr) DELETE samplesPtr;
			samplesPtr = NULL;

			//确定样本维度
			if ((this->sampleCount = samples.size()) <= 0) return false;
			this->sampleDimension = (*samples.begin())->Size(); //确定样本维度

			//创建样本数组
			this->samplesPtr = NEW("ksvdMatrix::CreateSamplesMatrix: Create Samples Matrix")
				Matrix<T>(this->sampleDimension, this->sampleCount);

			if (!samplesPtr) return false;

			//将样本压入样本矩阵
			Matrix<T> & samplesMatrix = *samplesPtr;

			auto i = samples.begin();
			auto end = samples.end(); //list<Vector<T> *>::iterator

			//枚举所有的样本
			for (int colIndex = 0; i != end; ++i, ++colIndex)
				samplesMatrix.SetCol(colIndex, **i);

			return true;
		}

		//创建码本数组
		bool CreateDictionaryMatrix(int codeBookCount)
		{
			this->codeBookCount = codeBookCount;				//确定码本个数

			//删除原有码本
			if (this->dictionaryPtr) DELETE this->dictionaryPtr;
			this->dictionaryPtr = NULL;

			//确定样本个数
			if ((this->sampleCount = samples.size()) <= 0) return false;
			this->sampleDimension = (*samples.begin())->Size();

			//创建新码本
			this->dictionaryPtr = NEW("ksvdMatrix::CreateDictionaryMatrix: Create Dictionary Matrix")
				Matrix<T>(this->sampleDimension, this->codeBookCount);

			if (dictionaryPtr) return true;
			return false;
		}

		bool CreateCoefficientMatrix(int codeBookCount)
		{
			this->codeBookCount = codeBookCount;				//确定码本个数

			//删除原有系数矩阵
			if (this->coefficientPtr) DELETE coefficientPtr;
			this->coefficientPtr = NULL;

			//确定样本个数
			if ((this->sampleCount = samples.size()) <= 0) return false;

			//创建系数矩阵
			this->coefficientPtr = NEW("ksvdMatrix::CreateCoefficientMatrix: Create Coefficient Matrix")
				Matrix<T>(this->codeBookCount, this->sampleCount);

			if (this->coefficientPtr) return true;
			return false;
		}

	public:
		//ksvd算法建立码本
		//codeBookCount: 输出码本的数量
		bool Run(int codeBookCount, int Sparsity, double eps, double Lemda,
			string saveFolderPath, int imgWidth, int imgHeight)
		{
			int beginI = 0;
			cout << "EPS: " << eps << endl;
			if (!this->CreateSamplesMatrix()) return false;					 //创建样本矩阵
			if (!this->CreateDictionaryMatrix(codeBookCount)) return false;	 //创建码本矩阵
			if (!this->CreateCoefficientMatrix(codeBookCount)) return false; //创建系数矩阵

			int n = this->sampleDimension;
			int N = this->sampleCount;
			int K = this->codeBookCount;
			cout << "Sample Count: " << sampleCount << endl;

			Matrix<T> & D = *this->dictionaryPtr;
			Matrix<T> & Y = *this->samplesPtr;
			Matrix<T> & X = *this->coefficientPtr;

			//将D初始化为样本的随机组合
			ifstream in(saveFolderPath + "\\Dictionary.txt");
			if (in.bad() || in.fail())
			{
				if (false)
				{
					D = RandomMatrix<T>(n, K, 0, 1);
					for (int i = 0; i < K; ++i)
					{
						Vector<T> v = D.GetCol(i);
						v.Normalization();
						D.SetCol(i, v);
					}
				}
				//D = RandomMatrix<T>(n, K, 0, 1);
				else
				{

					cout << "KMeas..." << endl;
					ca.kMeans(codeBookCount, nullptr, eps / 1E8, saveFolderPath, imgWidth, imgHeight);
					cout << "init Dictionary" << endl;
					for (int i = 0; i < K; ++i)
					{
						Vector<T> v = ca.Centers()[i];
						v.Normalization();
						D.SetCol(i, v);
					}
				}
			}
			else
			{
				D.ReadFromTextFile(in);
				in.close();
			}

			SaveDictionary(saveFolderPath, imgWidth, imgHeight);

			//D = RandomMatrix<T>(n, K, -1, 1);
			////cout << D;

			//将所有的初始码本规一化
			for (int i = 0; i < K; ++i)
			{
				Vector<T> column = D.GetCol(i);
				D.SetCol(i, column / column.Norm2());
			}
			//cout << "Dictionary: \n" << D << endl;
			Matrix<T> minD = D;
			double e = 0;
			int minId = 0;
			double minE = 0;
			eps = eps * D.AbsMax();
			for (int id = 0;; ++id)
			{

				//计算新的码本，更新系数矩阵X
				//cout << "OMP..." << endl;
				OMP<T>(D, Y, Sparsity, X, eps);
				////cout << "D\n" << D << endl;
				////cout << "Y\n" << Y << endl;
				////cout << "X\n" << X << endl;

				//计算差距
				Matrix<T> E = Y - D * X;
				////cout << "Err\n" << E << endl;

				double v1 = E.Norm2() / E.Col();
				double v2 = X.Norm0() * 1.0 / sampleCount;
				double checkValue = v1 + v2 * Lemda;
				cout << "EValue " << checkValue << " = " <<
					v1 << " + " << v2 << " * " << Lemda;
				if (id == 0)
				{
					minD = D;
					minE = checkValue;
					minId = 0;
					e = checkValue;
					cout << " min " << endl;
				}
				else
				{
					if (checkValue < minE)
					{
						minD = D;
						minId = id;
						minE = checkValue;
						cout << " min " << endl;
						if (checkValue < e && Abs<T>(e - checkValue) < eps / 1E2) break;
					}
					else cout << endl;

					e = checkValue;
				}

				//更新字典
				//int beginI = rand(K);
				//int IEnd = 1 + beginI;
				//for (int ii = beginI; ii < IEnd; ++ii)
				for (int i = 0; i < K; ++i)
				{
					//beginI++;
					//beginI = beginI % K;
					//int i = ii % K;

					//获得字典的第i行中的非0元的个数，以及索引
					Vector<T> xti = X.GetRow(i).Transpose();
					int wc = xti.Norm0(); //非0元个数
					bool empty = false;
					if (wc <= 0)
					{
						wc = 1;
						empty = true;
					}

					Vector<int> w(wc);
					int jj = 0;
					Matrix<real> Yj(n, wc);
					Matrix<real> Xj(K, wc);
					for (int j = 0; j < xti.Size(); ++j)
					{
						if (Abs(xti[j]) > epsilon || empty)
						{
							Yj.SetCol(jj, Y.GetCol(j));
							Xj.SetCol(jj, X.GetCol(j));
							w[jj++] = j;
							if (empty) break;
						}
					}



					D.SetCol(i, 0);				//清除字典的第i列
					X.SetRow(i, 0);				//清除系数的第i行
					Matrix<T> Ei = Y - D * X;
					Matrix<real> Eij(n, wc);
					for (int jj = 0; jj < wc; ++jj)
					{
						Eij.SetCol(jj, Ei.GetCol(w[jj]));
					}



					////cout << "Ei\n" << Ei << endl;

					//对Ei进行SVD分解
					Matrix<T> U, S, V;
					cout << "SVD ..." << i << '\r';
					Eij.SVD(U, S, V);
					T s0 = S[0][0];
					////cout << "S\n" << S << endl;
					D.SetCol(i, U.GetCol(0));
					Vector<T> nx = V.GetCol(0) * s0;
					VectorH<T> nix(N, 0);
					for (int jj = 0; jj < wc; ++jj)
					{
						nix[w[jj]] = nx[jj];
					}
					X.SetRow(i, nix);


					//SaveDictionary(saveFolderPath, imgWidth, imgHeight, i);

					//OMP<T>(D, Y, Sparsity, X, eps);
				}
				cout << endl;
				//cout << "D\n" << D << endl;

				SaveDictionary(saveFolderPath, imgWidth, imgHeight);
			}

			return true;
		}

		//获得码本向量
		Vector<T> GetCodebook(int codebookIndex)
		{
			if (!dictionaryPtr) return Vector<T>(1);
			return this->dictionaryPtr->GetCol(codebookIndex);
		}

	public:
		void SaveDictionary(string codePath, int imgWidth, int imgHeight, int index = -1)
		{
			Matrix<T> & DF = this->Dictionary();
			//码本写入文件
			if (index < 0)
			{
				string fn = codePath + "\\Dictionary.txt";
				ofstream sw(fn);
				sw << DF.Row() << endl;
				sw << DF.Col() << endl;
				sw << DF;
				sw.close();


				real min = DF.Min();
				real max = DF.Max();

				real k = Min(fabs(min), fabs(max));
				min = -k;
				max = k;

				for (int i = 0; i < this->codeBookCount; ++i)
				{
					Vector<real>  & centerVector = DF.GetCol(i);

					string fn = codePath + "\\" + i + ".bmp";
					xmv::Bitmap<real> bmp(centerVector.Buffer(), imgWidth, imgHeight);

					real min = bmp.Min();
					real max = bmp.Max();
					real k = Max(fabs(min), fabs(max));
					min = -k;
					max = k;

					bmp.Save(fn, min, max);

					//cout << "\rWriting ... " << i << "/" << centerCount;
				}
			}
			else
			{
				Vector<real>  & centerVector = DF.GetCol(index);
				real min = centerVector.Min();
				real max = centerVector.Max();
				real k = Max(fabs(min), fabs(max));
				min = -k;
				max = k;
				string fn = codePath + "\\" + index + ".bmp";
				xmv::Bitmap<real> bmp(centerVector.Buffer(), imgWidth, imgHeight);
				bmp.Save(fn, min, max);
			}
		}
	};
}