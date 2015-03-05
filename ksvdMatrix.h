#pragma once

#include <list>
#include <cstdlib>
#include "vector.h"
#include "matrix.h"
#include "common.h"
#include "memory.h"



namespace xmv
{
	///KSVD算法类
	template <typename T>
	class KSVDMatrix
	{
	private:
		list<Vector<T> *> samples;		//样本

		Matrix<T> * dictionaryPtr;		//字典矩阵：n x K，n为样本维度，K为码本数量（即最后要确定的码本）
		Matrix<T> * samplesPtr;			//样本矩阵：n x N, n为样本维度，N为样本数量
		Matrix<T> * coefficientPtr;		//系数矩阵：K x N, K为码本数量，N为样本数量
		int codeBookCount;				//K，码本数量
		int sampleCount;				//N，样本数量
		int sampleDimension;			//n, 样本维度

	public:
		KSVDMatrix()				//构造函数
			: dictionaryPtr(NULL), samplesPtr(NULL), coefficientPtr(NULL), 
			codeBookCount(0), sampleCount(0), sampleDimension(0)
		{
		}

		//析构函数
		~KSVDMatrix()					
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
		}

		//通过一个数组来添加一个样本。
		void AddSample(T * sampleArray, int arraySize)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(T * sampleArray, int arraySize)") 
				Vector<T>(sampleArray, arraySize));
		}

		//通过一个向量添加一个样本
		template<typename ST>
		void AddSample(Vector<ST> & sample)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(Vector<ST> & sample)") Vector<T>(sample));
		}

		//通过一个数组来添加一个样本。
		template<typename ST>
		void AddSample(ST * sampleArray, int arraySize)
		{
			this->samples.push_back(NEW("ClusterArrayer.AddSample(ST * sampleArray, int arraySize)") 
				Vector<T>(sampleArray, arraySize));
		}

		//删除所有的样本
		void ClearSamples()
		{
			list<Vector<T> *>::iterator i = samples.begin();
			list<Vector<T> *>::iterator end = samples.end();

			for(; i != end; ++i)
				DELETE(*i);

			samples.clear();
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
			for(int colIndex = 0; i != end; ++i, ++colIndex)
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
			this->dictionaryPtr = NEW("ksvdMatrix::CreateDictionaryMatrix: Create Dictionary Matrix") Matrix<T>(this->sampleDimension, this->codeBookCount);

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
			this->coefficientPtr = NEW("ksvdMatrix::CreateCoefficientMatrix: Create Coefficient Matrix") Matrix<T>(this->codeBookCount, this->sampleCount); 

			if (this->coefficientPtr) return true;
			return false;
		}

	private:
		//检测
		double Check()
		{
			if (this->coefficientPtr && this->samplesPtr && this->dictionaryPtr)
			{
				Matrix<T> checkMatrix = *this->samplesPtr - *this->dictionaryPtr * *this->coefficientPtr;
				return checkMatrix.SquareOfModulo();
			}
			return -1;
		}

	public:
		//ksvd算法建立码本
		//codeBookCount: 输出码本的数量
		bool KSVD(int codeBookCount, double eps)
		{
			if (!this->CreateSamplesMatrix()) return false;					 //创建样本矩阵
			if (!this->CreateDictionaryMatrix(codeBookCount)) return false;	 //创建码本矩阵
			if (!this->CreateCoefficientMatrix(codeBookCount)) return false; //创建系数矩阵

			for(double e = 0;;)
			{
				//计算新的码本


				//计算码本距离
				double checkValue = Check();
				if (Abs<T>(checkValue - e) < eps) break;
				e = checkValue;
			}

			return true;
		}

		//获得码本向量
		Vector<T> GetCodebook(int codebookIndex)
		{
			if (! dictionaryPtr) return Vector<T>(1);
			return this->dictionaryPtr->GetCol(codebookIndex);


		}
	};
}