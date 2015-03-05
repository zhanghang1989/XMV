#pragma once

#include <list>
#include <cstdlib>
#include "vector.h"
#include "common.h"
#include "memory.h"



namespace xmv
{
	template <typename T>
	class ClusterAnalyzer
	{
	private:
		Vector<T> * centers;
		int centerCount;
		list<Vector<T> *> samples;

	public:
		ClusterAnalyzer()
		{
			centers = 0;
			centerCount = 0;
		}

		~ClusterAnalyzer()
		{
			//DELETE samples
			ClearSamples();

			//DELETE centers
			if (centers) DELETEARR(centers);
		}

		int CenterCount() { return this->centerCount; }
		Vector<T> * Centers() { return this->centers; }



		void AddSample(Vector<T> & sample)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(Vector<T> & sample)") Vector<T>(sample));
		}

		void AddSample(T * sample, int size)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(T * sample, int size)") Vector<T>(sample, size));
		}

		template<typename ST>
		void AddSample(Vector<ST> & sample)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(Vector<ST> & sample)") Vector<T>(sample));
		}

		template<typename ST>
		void AddSample(ST * sample, int size)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(ST * sample, int size)") Vector<T>(sample, size));
		}

		void ClearSamples()
		{
			list<Vector<T> *>::iterator i = samples.begin();
			list<Vector<T> *>::iterator end = samples.end();

			for (; i != end; ++i)
				DELETE(*i);

			samples.clear();
		}

	private:
		//将样本压入数组buff。
		void FillSamplesBuffer(Vector<T> ** buff, int size)
		{
			list<Vector<T> *>::iterator i = samples.begin();
			list<Vector<T> *>::iterator end = samples.end();

			for (; i != end; ++i)
				*buff++ = *i;
		}

	private:
		void InitRandCenters(int centerCount, Vector<T> ** samplesPtr, int samplesCount)
		{
			if (centers) DELETEARR(centers);
			centers = NEW("ClusterAnalyzer.InitRandCenters: Centers") Vector<T>[centerCount];
			this->centerCount = centerCount;

			//准备随机数组
			int * temp = NEW("ClusterAnalyzer.InitRandCenters: Temp Rand Array") int[samplesCount];
			int * t = temp;
			for (int i = 0; i < samplesCount; ++i)
				*(t++) = i;

			//对随机数组前centerCount进行随机交换
			int * tEnd = temp + centerCount;
			int * t2, ts;
			for (t = temp; t < tEnd; ++t)
			{
				int k = (int)((long long)std::rand() * samplesCount / RAND_MAX);
				t2 = temp + k;
				ts = *t2;
				*t2 = *t;
				*t = ts;
			}

			//取前centerCount个作为初始聚类中心
			Vector<T>* sPtr;
			enumArray2(temp, centers, 0, centerCount, t, sPtr)
			{
				*sPtr = *samplesPtr[*t];
			} enumArrayEnd;

			DELETEARR temp;
		}

	private:
		Vector<T> *** CreateGroup(int groupCount, int groupCapacity)
		{
			Vector<T> *** group = NEW("ClusterAnalyzer.CreateGroup: Create Groups Array") Vector<T> **[groupCount];
			Vector<T> *** gPtr;
			enumArray(group, 0, groupCount, gPtr)
			{
				*gPtr = NEW("ClusterAnalyzer.CreateGroup: Create each group") Vector<T> *[groupCapacity];
				memset(*gPtr, 0, groupCapacity * sizeof(Vector<T> *));
			} enumArrayEnd;
			return group;
		}

	private:
		void ReleaseGroup(Vector<T>*** group, int groupCount)
		{
			Vector<T>*** gPtr;
			enumArray(group, 0, groupCount, gPtr)
			{
				DELETEARR(*gPtr);
			} enumArrayEnd;
			DELETEARR(group);
		}

	private:
		double KMeansGroup(Vector<T>*** group, int groupCount, Vector<T> ** samples, int sampleCount,
			int * groupIndexes)
		{
			Vector<T>*** groupEnds = NEW("ClusterAnalyzer.KMeansGroup: groupEnds") Vector<T>**[groupCount];
			memcpy(groupEnds, group, groupCount * sizeof(Vector<T>**));

			Vector<T>** sPtr;
			//为每一个样本计算距离
			double sum = 0;
			int sampleIndex = 0;
			enumArray(samples, 0, sampleCount, sPtr)
			{
				Vector<T> * cPtr, *nearestCenter = this->centers;
				Vector<T> & sample = **sPtr;
				double minDistance = this->centers->SquareOfDistance(sample);
				double distance = 0;
				//计算每一个样本和每一个中心的距离
				enumArray(this->centers, 1, groupCount, cPtr)
				{
					distance = Abs<T>((*cPtr - sample).SquareOfNorm2());
					if (distance < minDistance)
					{
						minDistance = distance;
						nearestCenter = cPtr;
					}
					////cout << distance << endl;
				} enumArrayEnd;

				sum += minDistance;
				//if (! (minDistance >= 0))
				//{
				//	//cout << sample << endl;
				//	//cout << (sPtr-samples) << " " << minDistance << "\t" << sum << endl;
				//	system("pause");
				//}
				//根据最近距离确定分组
				int groupID = nearestCenter - this->centers;
				*(groupEnds[groupID]++) = *sPtr;
				if (groupIndexes) groupIndexes[sampleIndex] = groupID;
				////cout << sum << endl;
				++sampleIndex;
			}enumArrayEnd;

			//将group末尾封闭
			Vector<T> *** gPtr;
			enumArray(groupEnds, 0, groupCount, gPtr)
			{
				**gPtr = NULL;
			} enumArrayEnd;

			DELETEARR groupEnds;

			//cout << sum << "\t";
			return sum;
		}

		double KMeansGroupCenter(Vector<T>*** group, int groupCount)
		{
			Vector<T>*** g;
			Vector<T>* c;
			Vector<T> sum;
			double distanceSum = 0;
			enumArray2(this->centers, group, 0, groupCount, c, g)
			{
				Vector<T> ** g2 = *g;
				if (*g2)
				{
					T * vPtr, *sumPtr;
					sum = **(g2++);
					int groupSize = 1;
					for (; *g2; ++g2, ++groupSize)
					{
						enumArray2((*g2)->Buffer(), sum.Buffer(), 0, (*g2)->Size(), vPtr, sumPtr)
							*sumPtr += *vPtr;
						enumArrayEnd;
					}
					sum /= groupSize;
					distanceSum += c->SquareOfDistance(sum);
					*c = sum;
				}
			} enumArrayEnd;

			return distanceSum;
		}

	public:
		Vector<T> * kMeans(int centerCount, int * groupIndexes = 0, double eps = EPS, string saveFolderPath = "", int imgWidth = 0, int imgHeight = 0)
		{
			//获得样本数组
			int sampleCount = this->samples.size();
			Vector<T> ** sampleBuff = NEW("ClusterAnalyzer.kMeans: Create Sample Buffer") Vector<T> *[sampleCount];
			this->FillSamplesBuffer(sampleBuff, sampleCount);

			//分组内存
			Vector<T> *** groups = this->CreateGroup(centerCount, sampleCount + 1);

			//随机选择centerCount个样本作为初始码本
			this->InitRandCenters(centerCount, sampleBuff, sampleCount);

			for (double d = eps * 2; d > eps;)
			{
				//根据现在的码本进行重新分组
				this->KMeansGroup(groups, centerCount, sampleBuff, sampleCount, groupIndexes);

				//计算每个分组的中心，并计算前后总体距离差
				d = this->KMeansGroupCenter(groups, centerCount);
				cout << d << endl;

				if (saveFolderPath != "")
					SaveDictionary(saveFolderPath, imgWidth, imgHeight);
			}

			//释放内存
			this->ReleaseGroup(groups, centerCount);
			DELETEARR(sampleBuff);

			return centers;
		}

	public:
		void SaveDictionary(string codePath, int imgWidth, int imgHeight)
		{
			//码本写入文件
//			string fn = codePath + "\\Dictionary.txt";
////			Matrix<T> DF = this->Dictionary();
//			ofstream sw(fn);
//			sw << DF.Row() << endl;
//			sw << DF.Col() << endl;
//			sw << DF;
//			sw.close();
//
//			real min = DF.Min();
//			real max = DF.Max();

			Vector<real> centerVector;

			for (int i = 0; i < this->centerCount; ++i)
			{
				centerVector = this->centers[i];

				string fn = codePath + "\\" + i + ".bmp";
				xmv::Bitmap<real> bmp(centerVector.Buffer(), imgWidth, imgHeight);
				bmp.Save(fn, bmp.Min(), bmp.Max());

				//cout << "\rWriting ... " << i << "/" << centerCount;
			}
		}

	private:
		void CreateAPSMatrix(Matrix<T> & s, Matrix<T> & ss)
		{
			//获得样本数组个数
			int sampleCount = this->samples.size();

			s.ReSize(sampleCount, sampleCount);
			int n = (*samples.begin())->Size();
			ss.ReSize(n, sampleCount);

			auto end = this->samples.end();

			T * p = s.Buffer();
			int id = 0;
			for (auto i = this->samples.begin(); i != end; ++i, ++id)
			{
				Vector<T> & si = **i;
				ss.SetCol(id, si);
				for (auto j = this->samples.begin(); j != end; ++j)
				{
					Vector<T> & sj = **j;
					*p++ = -(si - sj).Norm2();
				}
			}

			for (int i = 0; i < sampleCount; ++i)
			{
				Vector<T> k = s.GetCol(i);
				xmv::Quicksort(k.Buffer(), k.Size());
				s[i][i] = k[k.Size() / 2];
			}
		}


	public:
		Matrix<T> AP(int centerCount, int * groupIndexes = 0, double eps = EPS)
		{
			list<Vector<T>> centerList;
			//获得样本数组个数
			int sampleCount = this->samples.size();
			Vector<bool> K(sampleCount, true);

			//创建样本相似度矩阵S
			Matrix<T> S, mSamples;
			CreateAPSMatrix(S, mSamples);

			Matrix<T> R(sampleCount, sampleCount, 0);
			Matrix<T> A(sampleCount, sampleCount, 0);

			for (int id = 0; id < 100; ++id)
			{
				Matrix<T> newR(sampleCount, sampleCount, 0);
				Matrix<T> newA(sampleCount, sampleCount, 0);

				//R(k,k)=P(k)-max{A(k,j)+S(k,j)} (j {1,2,……,N,但j≠k})
				for (int k = 0; k < sampleCount; ++k)
				{
					VectorH<T> v = A.GetRow(k) + S.GetRow(k);
					v[k] = v.Min();
					newR[k][k] = S[k][k] - v.Max();
				}

				cout << "step 1 \n";

				//A(i,k)=min{0, R(k,k)+  Sigma(j)(max(0, R(j, k))) } {j = 1,2,……,N,但j≠i且j≠k}) 
				for (int k = 0; k < sampleCount; ++k)
				{
					Vector<double> v = R.GetCol(k);
					for (int i = 0; i < sampleCount; ++i)
					{
						double sum = 0;
						for (int j = 0; j < v.Size(); ++j)
						if (j != i) sum += max(0, v[j]);
						double aik = newR[k][k] + sum;
						newA[i][k] = min(0, aik);
					}
				}

				cout << "step 2 \n";
				//R(i,k)=S(i,k)- max{A(i,j)+S(i,j)}(j {1,2,……,N,但j≠k})
				for (int k = 0; k < sampleCount; ++k)
				{
					for (int i = 0; i < sampleCount; ++i)
					{
						VectorH<T> v = newA.GetRow(i) + S.GetRow(i);
						v[k] = v.Min();
						newR[i][k] = v.Max() + S[i][k];
					}
				}

				cout << "step 3 \n";

				double lam = 0.75;
				R = (1 - lam) * R + lam * newR;
				A = (1 - lam) * A + lam * newA;
			}

			for (int k = 0; k < sampleCount; ++k)
			{
				if (R[k][k] + A[k][k] > 0)
				{
					centerList.push_back(mSamples.GetCol(k));
					cout << mSamples.GetCol(k) << endl;
				}
			}

			Matrix<T> mOut(mSamples.Row(), centerList.size());
			auto end = centerList.end();
			int id = 0;
			for (auto i = centerList.begin(); i != end; ++i, ++id)
				mOut.SetCol(id, *i);
			return mOut;
		}
	};

}