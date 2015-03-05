#pragma once
#include <cstring>
#include <cstdarg>
#include <iostream>
#include <cstdlib>
#include "matrix.h"
#include "memory.h"
#include "array.h"

using namespace std;

namespace xmv
{
	template<typename T>
	class VectorH : public xmv::Array<T>
	{

	public:
		VectorH();
		VectorH(int size);
		VectorH(int size, T value);
		VectorH(VectorH & other);
		VectorH(T * arr, int size);
		template <typename T2> VectorH(T2 * arr, int size);
		template <typename T2> VectorH(VectorH<T2> & other);

		~VectorH(void);

		Vector<T> Transpose(void);
		Vector<T> ConjugateTranspose(void);

		//从向量中截取一部分成为新的向量
		VectorH<T> Cut(int start, int len)
		{
			VectorH<T> v(len);
			memcpy(v.dataPtr, this->dataPtr + start, len * sizeof(T));
			return v;
		}
		//从向量中截取一部分成为新的向量
		VectorH<T> Cut(int start)
		{
			int len = this->size - start;
			VectorH<T> v(len);
			memcpy(v.dataPtr, this->dataPtr + start, len * sizeof(T));
			return v;
		}


	public:

		inline VectorH<T> & operator = (VectorH<T> & v);
		inline VectorH<T> & operator = (T value);
		inline VectorH<T> & operator += (VectorH<T> & v);
		inline VectorH<T> & operator -= (VectorH<T> & v);
		inline VectorH<T> & operator *= (T value);
		inline VectorH<T> & operator /= (T value);
		inline T & operator[] (int index);
		inline Vector<T> operator ~();

		template<typename T2> VectorH<T> & operator = (VectorH<T2> & v);


		template<typename T> friend VectorH<T> operator + (VectorH<T> & v1, VectorH<T> & v2);
		template<typename T> friend VectorH<T> operator - (VectorH<T> & v1, VectorH<T> & v2);
		template<typename T> friend VectorH<T> operator * (VectorH<T> & v1, Matrix<T> & v2);
		template<typename T> friend VectorH<T> operator + (VectorH<T> & v);
		template<typename T> friend VectorH<T> operator - (VectorH<T> & v);
		template<typename T> friend VectorH<T> operator * (T value, VectorH<T> & vectorH);
		template<typename T> friend VectorH<T> operator * (VectorH<T> & vectorH, T value);
		template<typename T> friend VectorH<T> operator / (VectorH<T> & vectorH, T value);
		template<typename T> friend T operator * (VectorH<T> & v1, VectorH<T> & v2);

		template<typename T> friend ostream & operator << (ostream & out, VectorH<T> & vectorH);
		template<typename T> friend istream & operator >> (istream & in, VectorH<T> & vectorH);


	};

	template<typename T>
	VectorH<T>::VectorH() : Array() 
	{
	}

	template<typename T>
	VectorH<T>::VectorH(int size) : Array(size)
	{
	}

	template<typename T>
	VectorH<T>::VectorH(VectorH<T> & other) : Array(other.dataPtr, other.size)
	{
	}

	template<typename T>
	template<typename T2>
	VectorH<T>::VectorH(VectorH<T2> & other) : Array(other.Buffer(), other.Size())
	{
	}

	template<typename T>
	VectorH<T>::VectorH(int size, T value) : Array(size, value)
	{
	}

	template<typename T>
	VectorH<T>::VectorH(T * arr, int size) : Array(arr, size)
	{
	}

	template<typename T>
	template<typename T2>
	VectorH<T>::	VectorH(T2 * arr, int size) : Array(arr, size)
	{
	}

	template<typename T> Vector<T> VectorH<T>::Transpose()
	{
		Vector<T> m(size);
		memcpy(m.Buffer(), this->dataPtr, this->memorySize);
		return m;
	}

	template<typename T> Vector<T> VectorH<T>::ConjugateTranspose()
	{
		Vector<T> m(size);
		T * ptr1, * ptr2;
		enumArray2(this->dataPtr, m.Buffer(), 0, size, ptr1, ptr2)
			*ptr2 = Conjugate(*ptr1);
		enumArrayEnd;
		return m;
	}

	template<typename T> VectorH<T>::~VectorH(void)
	{
	}

	//extern DigtalFormat ____currentDigtalFormat;

	template<typename T> ostream & operator << (ostream & out, VectorH<T> & vectorH)
	{
		T * ptr = vectorH.dataPtr;
		int size = vectorH.size;

		out << ____currentDigtalFormat;

		for(int y = 0; y < size; ++y)
		{
			out << setw(____currentDigtalFormat.space) << *ptr++ << ' ';
		}

		return out;
	}

	template<typename T> istream & operator >> (istream & in, VectorH<T> & vectorH)
	{
		T * ptr = vectorH.dataPtr;
		int size = vectorH.size;

		for(int y = 0; !in.eof() && y < size; ++y)
		{
			in >> *ptr++;
		}

		return in;
	}

	template<typename T> inline Vector<T> VectorH<T>::operator ~ ()
	{
		return this->ConjugateTranspose();
	}

	template<typename T> VectorH<T> & VectorH<T>::operator = (VectorH<T> & v)
	{
		this->ReSize(v.size);

		memcpy(this->dataPtr, v.dataPtr, this->memorySize);

		return *this;
	}

	template<typename T> template<typename T2> VectorH<T> & VectorH<T>::operator = (VectorH<T2> & v)
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

	template<typename T> VectorH<T> & VectorH<T>::operator = (T value)
	{
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr)
			*ptr = value;
		return *this;
	}

	template<typename T> VectorH<T> & VectorH<T>::operator += (VectorH<T> & v)
	{
		T * sourcePtr = v.Buffer();
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr, ++sourcePtr)
			*ptr += *sourcePtr;
		return *this;
	}

	template<typename T> VectorH<T> & VectorH<T>::operator -= (VectorH<T> & v)
	{
		T * sourcePtr = v.Buffer();
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr, ++sourcePtr)
			*ptr -= *sourcePtr;
		return *this;
	}

	template<typename T> VectorH<T> & VectorH<T>::operator *= (T value)
	{
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr)
			*ptr *= value;
		return *this;
	}

	template<typename T> VectorH<T> & VectorH<T>::operator /= (T value)
	{
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr)
			*ptr /= value;
		return *this;
	}

	template<typename T> VectorH<T> operator + (VectorH<T> & v1, VectorH<T> & v2)
	{
		VectorH<T> v = v1;
		v += v2;
		return v;
	}

	template<typename T> VectorH<T> operator - (VectorH<T> & v1, VectorH<T> & v2)
	{
		VectorH<T> v = v1;
		v -= v2;
		return v;
	}

	template<typename T> VectorH<T> operator + (VectorH<T> & vectorH)
	{
		VectorH<T> v = VectorH;
		return v;
	}

	template<typename T> VectorH<T> operator - (VectorH<T> & vectorH)
	{
		VectorH<T> v(VectorH.size);
		T * sourcePtr = VectorH.dataPtr;
		for(T * ptr = v.dataPtr; ptr < v.dataEndPtr; ++ptr, ++sourcePtr)
			*ptr = -*sourcePtr;
		return v;
	}

	template<typename T> VectorH<T> operator * (VectorH<T> & m1, Matrix<T> & m2)
	{
		int c1 = m1.size;
		int r1 = 1;
		int c2 = m2.Col();
		int r2 = m2.Row();
		if (c1 != r2) throw;
		VectorH<T> m(c2, 0);
		int x, y;
		T * m1sizePtr = m1.dataPtr;
		T * ptr1;
		T * m2Ptr = m2.Buffer(), *m2ColPtr, *ptr2;
		T * ptr = m.Buffer();

		for(x = 0; x < c2; ++x, ++ptr) //循环m的行
		{
			ptr1 = m1.Buffer();
			ptr2 = m2.Buffer() + x;
			for(y = 0; y < r2; ++y, ++ptr1, ptr2 += c2) //循环m的列
			{
				*ptr += *ptr1 * *ptr2;
			}
		}

		return m;
	}

	template<typename T> VectorH<T> operator * (T val, VectorH<T> & vectorH)
	{
		VectorH<T> m = VectorH;
		m *= val;

		return m;
	}

	template<typename T> VectorH<T> operator * (VectorH<T> & vectorH, T val)
	{
		VectorH<T> m = VectorH;
		m *= val;

		return m;
	}

	template<typename T> VectorH<T> operator / (VectorH<T> & vectorH, T val)
	{
		VectorH<T> m = VectorH;
		m /= val;

		return m;
	}

	template<typename T> T & VectorH<T>::operator[] (int index)
	{
		return this->dataPtr[index];
	}

	template<typename T> T operator * (VectorH<T> & v1, VectorH<T> & v2)
	{
		if (v1.Size() != v2.Size()) throw;

		T * p1, * p2;
		T sum = 0;
		enumArray2(v1.Buffer(), v2.Buffer(), 0, v1.Size(), p1, p2)
			sum += *p1 * *p2;
		enumArrayEnd;

		return sum;
	}
};