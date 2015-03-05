#pragma once
#include <cstring>
#include <cstdarg>
#include <iostream>
#include <cstdlib>

#include "memory_.h"
#include "convert.h"
#include "common.h"
#include "io.h"
#include <fstream>
#include "rgb.h"

using namespace std;

#pragma once

extern double epsilon;

namespace xmv
{
#define EnumArr(type, pointer) for(T * pointer = dataPtr; pointer < dataEndPtr; ++pointer)


	template<typename T>
	class Array
	{
	protected:
		int size;								//size count of the Array
		int memorySize;							//memory size of the Array  = size * sizeof(T)
		T * dataPtr;							//the start fo Array array
		T * dataEndPtr;							//the end of Array array

	public:
		Array();
		Array(int size);
		Array(int size, T value);
		Array(Array<T> & other);
		Array(T * arr, int size);
		Array(BGR<T> * arr, int size);
		Array(int size, T* arr, int arrSize);
		template <typename T2> Array(Array<T2> & other);
		template <typename T2> Array(T2 * arr, int size);
		template <typename T2> Array(int size, T2 * arr, int arrSize);

		~Array(void);

		inline int Size();						// 返回元素的个数
		inline int MemorySize();				// 返回实际占用内存的字节数
		inline void ReSize(int size);			// 重新设定大小，清除所有的数据
		inline T * Buffer();					// 获取内部数组的首地址
		inline void SetValue(int index, T val);	// 设置第index个数据的值
		inline T GetValue(int index);			// 获得第index个数据
		T SquareOfNorm2();						// 2范式的平方
		T SquareOfDistance(Array<T> & another);	// 距离另一个数组的距离的平方
		T Distance(Array<T> & another);			// 计算两个数组的距离
		Array<T> Cut(int start, int len);		// 从向量中截取一部分成为新的向量
		Array<T> Cut(int start);				// 将数组剪切一部分，返回新的数组，从start开始，一直到结束
		T Norm0();								// 0范式，即求非0元素个数
		T Norm1();								// 1范式，即求所有元素的和
		T Norm2();								// 2范式，即求所有元素的平方和开根号
		T Norm(int p = 2);						// p范式
		T Sum();								// 所有元素的和
		void SetSmallValueToZero();				// 将绝对值小于eps的数值设置为0
		T AbsMin();								// 返回绝对值的最小值
		T Min();								// 返回最小值
		T AbsMax();								// 返回绝对值的最大值
		T Max();								// 返回最大值
		void Offset(T value);					// 所有数据偏移value
		T Average();							// 获得所有数据的平均值
		void AverageTo(T value = 0);			// 将数据平均化到某个数值
		T Variance();							// 返回方差
		void GaussianWhiten();					// 高斯白化，使得数组平均值为0，方差为1
		void AbsValues();						// 所有的元素取绝对值
		bool InfinityCheck(T maxValue = 1E50);	// 数据检查，看数组中是否有数据为无穷大
		int ClearInfinityValues(T maxValue = 1E50, T clearValue = 0);	// 数据检查，看数组中是否有数据为无穷大，如果有，则清空为clearvalue，返回清空的数据的个数
		bool IsZero();							// 返回数据是否全部为0
		void Clear();							// 将数据全部清0
		void Normalization();					// 正规化，每个元素平方和为1
		void SumAsOne();						// 是每个元素的和为1
		void SetMinTo(T value);					// 将小于value的元素全部设置为value
		void SetMaxTo(T value);					// 将大于value的元素全部设置为value
		bool Contains(T value, int begin = 0, int len = -1); //返回是否包含值为value的元素=-0987654321`
		Array<int> histgram(int binnum);

		virtual void WriteToTextFile(string filename, bool append = false);
		virtual void WriteToTextFile(ofstream & o);
		virtual void ReadFromTextFile(string filename);
		virtual void ReadFromTextFile(ifstream & i);

	protected: 
		inline void Init(int size);

	public:
		Array<T> & operator = (Array<T> & v);
		Array<T> & operator = (T value);
		T & operator[] (int index);
		template<typename T2> Array<T> & operator = (Array<T2> & v);

		template<typename T> friend ostream & operator << (ostream & out, Array<T> & vector);
		template<typename T> friend istream & operator >> (istream & in, Array<T> & vector);
	};

	template<typename T>
	inline int Array<T>::Size()
	{
		return size; 
	}

	template<typename T>
	inline int Array<T>::MemorySize()
	{ 
		return memorySize;
	}

	template<typename T>
	inline T * xmv::Array<T>::Buffer()
	{
		return dataPtr; 
	}

	template<typename T>
	inline void xmv::Array<T>::SetValue(int index, T val)
	{
		dataPtr[index] = val;
	}

	template<typename T>
	inline T Array<T>::GetValue(int index)
	{
		return dataPtr[index]; 
	}

	//高斯白化
	template<typename T> void Array<T>::GaussianWhiten()
	{
		AverageTo(0);
		T val = Sqrt<T>(this->Variance());
		if (val == 0) 
		{
			val = Sqrt<T>(this->Average());
			return;
		}

		EnumArr(T, ptr)
			*ptr /= val;
	}

	template<typename T>
	T Array<T>::SquareOfNorm2()
	{
		T sum = 0;
		EnumArr(T, t)
			sum += *t * *t;
		return sum;
	}

	template<typename T>
	T Array<T>::Norm2()
	{ 
		return Sqrt<T>(SquareOfNorm2()); 
	}

	template<typename T>
	T Array<T>::SquareOfDistance(Array<T> & another)
	{
		T sum = 0;
		T * p1, *p2, t;
		int s = xmv::Min(this->size, another.size);
		enumArray2(this->dataPtr, another.dataPtr, 0, s, p1, p2)
		{
			t = *p1 - *p2;
			sum += t * t;
		} enumArrayEnd;
		return sum;
	}


	template<typename T>
	T Array<T>::Distance(Array<T> & another)
	{
		return Sqrt<T>(SquareOfDistance(another));
	}

	template<typename T>
	Array<T> Array<T>::Cut(int start, int len)
	{
		Array<T> v(len);
		memcpy(v.dataPtr, this->dataPtr + start, len * sizeof(T));
		return v;
	}

	template<typename T>
	Array<T> Array<T>::Cut(int start)
	{
		int len = this->size - start;
		Array<T> v(len);
		memcpy(v.dataPtr, this->dataPtr + start, len * sizeof(T));
		return v;
	}


	template<typename T> T Array<T>::Average()					// 获得所有数据的平均值
	{
		T sum = 0;
		EnumArr(T, ptr)
			sum += *ptr;
		return sum / size;
	}

	template<typename T> T Array<T>::Variance()							// 返回方差
	{
		T avg = Average();
		T sum = 0;
		T temp;
		EnumArr(T, ptr)
		{
			temp = *ptr - avg;
			sum += temp * temp;
		}
		return sum / this->size;
	}

	template<typename T>
	void Array<T>::AverageTo(T value = 0)
	{
		T avg = this->Average();
		EnumArr(T, ptr)
		{
			*ptr -= (avg - value);
		}
	}

	template<typename T>
	T Array<T>::Norm0() // 0范式，即求非0元素个数
	{
		int count = 0;
		EnumArr(T, ptr)
			if (Abs<T>(*ptr) > epsilon) count++;
		return count;
	}

	template<typename T>
	T Array<T>::Norm1() // 1范式，即求所有元素的和
	{
		T sum = 0;
		EnumArr(T, ptr)
			sum += *ptr;
		return sum;
	}

	template<typename T>
	T Array<T>::Sum() // 即求所有元素的和
	{
		T sum = 0;
		EnumArr(T, ptr)
			sum += *ptr;
		return sum;
	}

	template<typename T>
	inline void Array<T>::Init(int size)
	{
		this->size = size;
		memorySize = size * sizeof(T);

		dataPtr = NEW("XMV.Array<T>.Init: Create Array Buffer") T[size];
		dataEndPtr = dataPtr + size;
	}

	template<typename T>
	Array<T>::Array()
	{ 
		Init(1); 
	}

	template<typename T>
	Array<T>::Array(int size)
	{
		Init(size);
	}

	template<typename T>
	Array<T>::Array(int size, T value)
	{
		Init(size);

		EnumArr(T, ptr)
			*ptr = value;
	}

	template<typename T>
	Array<T>::	Array(T * arr, int size)
	{
		Init(size);
		memcpy(this->dataPtr, arr, this->memorySize);
	}

	template<typename T>
	Array<T>::Array(int size, T* arr, int arrSize)
	{
		Init(size);
		memcpy(this->dataPtr, arr, min(this->memorySize, arrSize * sizeof(T)));
	}

	template<typename T>
	template<typename T2>
	Array<T>::	Array(T2 * arr, int size)
	{
		Init(size);
		T * p1;
		T2 * p2;
		enumArray2(this->dataPtr, arr, 0, size, p1, p2)
		{
			*p1 = *p2;
		}
		enumArrayEnd;
	}

	template<typename T>
	Array<T>::	Array(BGR<T> * arr, int size)
	{
		Init(size);
		T * p1;
		BGR<T> * p2;
		enumArray2(this->dataPtr, arr, 0, size, p1, p2)
		{
			*p1 = p2->Gray();
		}
		enumArrayEnd;
	}

	template<typename T>
	template<typename T2>
	Array<T>::	Array(int size, T2* arr, int arrSize)
	{
		Init(size);
		T * p1;
		T2 * p2;
		enumArray2(this->dataPtr, arr, 0, min(size, arrSize), p1, p2)
		{
			*p1 = *p2;
		}
		enumArrayEnd;
	}

	template<typename T>
	Array<T>::Array(Array<T> & other)
	{
		Init(other.size);
		::memcpy(this->dataPtr, other.dataPtr, memorySize);
	}

	template<typename T>
	template<typename T2>
	Array<T>::Array(Array<T2> & other)
	{
		Init(other.Size());
		T * ptr;
		T2 * ptr2;
		enumArray2(this->dataPtr, ohter.Buffer(), 0, size, ptr, ptr2)
			*ptr = (*ptr2);
		enumArrayEnd;
	}

	template<typename T>
	inline void Array<T>::ReSize(int size) 
	{
		if (this->size != size)
		{
			DELETEARR(dataPtr);
			Init(size);
		}
	}

	template<typename T> Array<T>::~Array(void)
	{
		DELETEARR(dataPtr);
	}

	extern DigtalFormat ____currentDigtalFormat;

	template<typename T> ostream & operator << (ostream & out, Array<T> & data)
	{
		out << ____currentDigtalFormat;

		T * ptr;
		for(T * ptr = data.dataPtr; ptr < data.dataEndPtr; ++ptr)
		{
			out << setw(____currentDigtalFormat.space) << *ptr << ' ';
		}

		return out;
	}

	template<typename T> istream & operator >> (istream & in, Array<T> & data)
	{
		T * ptr = data.dataPtr;
		int size = data.size;

		for(int y = 0; !in.eof() && y < size; ++y)
		{
			in >> *ptr++;
		}

		return in;
	}

	template<typename T> Array<T> & Array<T>::operator = (Array<T> & v)
	{
		this->ReSize(v.size);

		memcpy(this->dataPtr, v.dataPtr, this->memorySize);

		return *this;
	}

	template<typename T> template<typename T2> Array<T> & Array<T>::operator = (Array<T2> & v)
	{
		this->ReSize(v.Size());
		T * p1;
		T2 * p2;
		enumArray2(this->dataPtr, v.Buffer(), 0, this->size, p1, p2)
		{
			*p1 = (T)(*p2);
		}enumArrayEnd

			return *this;
	}

	template<typename T> Array<T> & Array<T>::operator = (T value)
	{
		EnumArr(T, ptr)
			*ptr = value;
		return *this;
	}

	template<typename T> T & Array<T>::operator[] (int index)
	{
		return this->dataPtr[index];
	}

	//将小于epsilon的值设置为0
	template<typename T> void Array<T>::SetSmallValueToZero() 
	{ 
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr)
			if (Abs<T>(*ptr) < epsilon) *ptr = 0;
		return *this;
	}

	template<typename T> T Array<T>::AbsMin()
	{
		T max;
		max = 0;
		EnumArr(T, ptr)
			if (Abs<T>(*ptr) < max) max = Abs<T>(*ptr);
		return max;
	}
	template<typename T> T Array<T>:: Min()
	{
		T max;
		max = 0;
		EnumArr(T, ptr)
			if (*ptr < max) max = *ptr;
		return max;
	}
	template<typename T> T Array<T>::AbsMax()
	{
		T max;
		max = 0;
		EnumArr(T, ptr)
			if (Abs<T>(*ptr) > max) max = Abs<T>(*ptr);
		return max;
	}

	template<typename T> T Array<T>::Max()
	{
		T max;
		max = 0;
		EnumArr(T, ptr)
			if (*ptr > max) max = *ptr;
		return max;
	}

	template<typename T> void Array<T>::Offset(T value)
	{
		T max = 0;
		EnumArr(T, ptr)
			* ptr += value;
	}

	template<typename T> T Array<T>::Norm(int p)	//范式
	{
		T sum = 0;
		EnumArr(T, ptr)
			sum += Pow<T>(*ptr, p);
		return Pow<T>(sum, 1.0 / p);
	}		

	template<typename T> void Array<T>::AbsValues()
	{
		EnumArr(T, ptr)
			*ptr = Abs<T>(*ptr);
	}

	// 数据检查，看数组中是否有数据为无穷大
	template<typename T> bool Array<T>::InfinityCheck(T maxValue)		// 数据检查，看数组中是否有数据为无穷大
	{
		EnumArr(T, ptr)
		{
			if (! (*ptr > -maxValue)) return false;
			if (! (*ptr < maxValue)) return false;
		}
		return true;
	}

	// 数据检查，看数组中是否有数据为无穷大，如果有，则清空为clearvalue，返回清空的数据的个数
	template<typename T> int Array<T>::ClearInfinityValues(T maxValue, T clearValue)
	{
		int count = 0;
		EnumArr(T, ptr)
		{
			if (! (*ptr > -maxValue)) { count++; *ptr = clearValue; }
			if (! (*ptr < maxValue)) { count++; *ptr = clearValue; }
		}
		return count;
	}


	// 数据检查，看数组中是否有数据为无穷大
	template<typename T> bool Array<T>::IsZero()		// 数据检查，看数组中是否有数据为无穷大
	{
		EnumArr(T, ptr)
		{
			if (*ptr > epsilon) return false;
			if (*ptr < -epsilon) return false;
		}
		return true;
	}

	//将数据清0
	template<typename T> void Array<T>::Clear()
	{
		memset(dataPtr, 0, memorySize);
	}

	// 正规化，每个元素平方和为1
	template<typename T> void Array<T>::Normalization()
	{
		T v = this->Norm2();
		if (Abs<T>(v) > epsilon)
			EnumArr(T, ptr)
			*ptr /= v;
	}

	// 使每个元素的和为1
	template<typename T> void Array<T>::SumAsOne()
	{
		T v = this->Norm1();
		if (Abs<T>(v) > epsilon)
			EnumArr(T, ptr)
			*ptr /= v;
	}

	// 将小于value的元素全部设置为value
	template<typename T> void Array<T>::SetMinTo(T value)
	{
		EnumArr(T, ptr)
			if (*ptr < value) *ptr = value;
	}
	// 将大于value的元素全部设置为value
	template<typename T> void Array<T>::SetMaxTo(T value)
	{			
		EnumArr(T, ptr)
			if (*ptr > value) *ptr = value;
	}


	template<typename T> void Array<T>::WriteToTextFile(ofstream & o)
	{
		o << this->size << endl;
		o << format(10, 12) << *this;
	}

	template<typename T> void Array<T>::WriteToTextFile(string filename, bool append = false)
	{
		ofstream o;

		if (append)
			o.open(filename, std::ios_base::ate);
		else
			o.open(filename);

		this->WriteToTextFile(o);

		o.close();
	}

	template<typename T> void Array<T>::ReadFromTextFile(string filename)
	{
		ifstream i;
		this->ReadFromTextFile(i);
		i.close();
	}

	template<typename T> void Array<T>::ReadFromTextFile(ifstream & i)
	{
		int s;
		i >> s;
		ReSize(s);

		i >> *this;
	}

	//返回是否包含值为value的元素
	template<typename T> bool Array<T>::Contains(T value, int begin = 0, int len = -1)
	{
		if (len < 0) len = size - begin;
		T * p;
		enumArray(this->dataPtr, begin, begin + len, p)
		{
			if (*p == value) 
				return true;
		}
		enumArrayEnd;

		return false;
	}

	// get the histgram of array
	template <typename T> Array<int> Array<T>::histgram(int binnum)
	{
		T step = (Max() - Min()) / binnum;
		Array<int> hist(binnum, 0);

		for (T* ptr = Buffer(); ptr < dataEndPtr; ptr++)
		{
			hist[(int) (*(ptr) / step) % binnum] ++;
		}
		return hist;
	}
}