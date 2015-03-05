#pragma once;

#include <stdlib.h>
#include <list>
using namespace std;

#include "array.h"
#include "vector.h"
#include "matrix.h"
#include "vectorH.h"



namespace xmv
{

	template<typename T> class SampleLoader
	{
	private:
		list<Vector<T>> samples;		//样本     
		int sampleSize;

	public:
		SampleLoader() : sampleSize(0) { }

		int SampleCount() { return samples.size(); }

		void AddSample(Vector<T> v) { sampleSize = v.Size(); samples.push_back(v); }

		void AddSample(T * arr, int size) { AddSample(Vector<T>(arr, size)); }

		//获得样本矩阵，其行为样本维度，列为样本个数，也就是说，每一列为一个样本
		Matrix<T> GetSampleMatrix()
		{
			Matrix<T> m(sampleSize, samples.size());

			auto end = samples.end();
			int id = 0;
			for(auto i = samples.begin(); i != end; ++i, ++id)
				m.SetCol(id, *i);

			return m;
		}
	};
	
}