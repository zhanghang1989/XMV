#pragma once
#include <cstring>
#include <cstdarg>
#include <iostream>
#include <cstdlib>
#include "memory.h"
#include "vectorH.h"
#include "array.h"
#include "io.h"

using namespace std;

namespace xmv
{

	template<typename T> class Vector : public xmv::Array<T>
	{
	public:
		Vector();
		Vector(int size);
		Vector(int size, T value);
		Vector(const Vector<T> & other);
		Vector(T * arr, int size);
		template <typename T2> Vector(T2 * arr, int size);
		template <typename T2> Vector(Vector<T2> & other);

		~Vector(void);

		inline VectorH<T> Transpose(void);				//返回该向量的转置
		inline VectorH<T> ConjugateTranspose(void);		//返回该向量的共轭转置
		inline Vector<T> Cut(int start, int len);		//从向量中截取一部分成为新的向量
		inline Vector<T> Cut(int start);				//从向量中截取一部分成为新的向量

	public:

		inline Vector<T> & operator = (Vector<T> & v);
		inline Vector<T> & operator = (T value);
		inline Vector<T> & operator += (Vector<T> & v);
		inline Vector<T> & operator -= (Vector<T> & v);
		inline Vector<T> & operator *= (T value);
		inline Vector<T> & operator /= (T value);
		inline T & operator[] (int index);
		inline VectorH<T> operator ~ ();

		template<typename T2> inline Vector<T> & operator = (Vector<T2> & v);

		template<typename T> friend Vector<T> operator + (Vector<T> & v1, Vector<T> & v2);
		template<typename T> friend Vector<T> operator - (Vector<T> & v1, Vector<T> & v2);
		template<typename T> friend Matrix<T> operator * (Vector<T> & v1, Matrix<T> & v2);
		template<typename T> friend Vector<T> operator * (Matrix<T> & v1, Vector<T> & v2);
		template<typename T> friend Vector<T> operator + (Vector<T> & v);
		template<typename T> friend Vector<T> operator - (Vector<T> & v);
		template<typename T> friend Vector<T> operator * (T value, Vector<T> & vector);
		template<typename T> friend Vector<T> operator * (Vector<T> & vector, T value);
		template<typename T> friend T operator * (Vector<T> & v1, Vector<T> & v2);
		template<typename T> friend T operator * (VectorH<T> & v1, Vector<T> & v2);
		template<typename T> friend Matrix<T> operator * (Vector<T> & v1, VectorH<T> & v2);
		template<typename T> friend Vector<T> operator / (Vector<T> & vector, T value);

		template<typename T> friend ostream & operator << (ostream & out, Vector<T> & vector);
		template<typename T> friend istream & operator >> (istream & in, Vector<T> & vector);
	};

	template<typename T>
	Vector<T>::Vector(const Vector<T> & other) : Array(other.dataPtr, other.size)
	{
	}

	template<typename T>
	template<typename T2>
	Vector<T>::Vector(Vector<T2> & other) : Array(other.Buffer(), other.Size())
	{
	}

	template<typename T>
	Vector<T>::Vector() : Array()
	{
	}

	template<typename T>
	Vector<T>::Vector(int size) : Array(size)
	{
	}

	template<typename T>
	Vector<T>::Vector(int size, T value) : Array(size, value)
	{
	}

	template<typename T>
	Vector<T>::Vector(T * arr, int size) : Array(arr, size)
	{
	}

	template<typename T>
	template<typename T2>
	Vector<T>::Vector(T2 * arr, int size) : Array(arr, size)
	{
	}

	template<typename T> Vector<T>::~Vector(void)
	{
	}

	template<typename T> VectorH<T> Vector<T>::Transpose()
	{
		VectorH<T> m(size);
		memcpy(m.Buffer(), this->dataPtr, this->memorySize);
		return m;
	}

	template<typename T> inline VectorH<T> Vector<T>::ConjugateTranspose()
	{
		VectorH<T> m(size);
		T * ptr1, * ptr2;
		enumArray2(this->dataPtr, m.Buffer(), 0, size, ptr1, ptr2)
			*ptr2 = Conjugate(*ptr1);
		enumArrayEnd;
		return m;
	}

	template<typename T> inline VectorH<T> Vector<T>::operator ~ ()
	{
		return this->ConjugateTranspose();
	}

	extern DigtalFormat ____currentDigtalFormat;

	template<typename T> ostream & operator << (ostream & out, Vector<T> & vector)
	{
		T * ptr = vector.dataPtr;
		int size = vector.size;

		out << ____currentDigtalFormat;


		for(int y = 0; y < size; ++y)
		{
			out << setw(____currentDigtalFormat.space) << *ptr++ << ' ';
		}

		return out;
	}

	template<typename T> istream & operator >> (istream & in, Vector<T> & vector)
	{
		T * ptr = vector.dataPtr;
		int size = vector.size;

		for(int y = 0; !in.eof() && y < size; ++y)
		{
			in >> *ptr++;
		}

		return in;
	}

	template<typename T> Vector<T> & Vector<T>::operator = (Vector<T> & v)
	{
		this->ReSize(v.size);

		memcpy(this->dataPtr, v.dataPtr, v.memorySize);

		return *this;
	}

	template<typename T> template<typename T2> Vector<T> & Vector<T>::operator = (Vector<T2> & v)
	{
		this->ReSize(v.Size());
		T * p1;
		T2 * p2;
		enumArray2(this->dataPtr, v.Buffer(), 0, this->size, p1, p2)
		{
			*p1 = (T)(*p2);
		}enumArrayEnd;

		return *this;
	}

	template<typename T> Vector<T> & Vector<T>::operator = (T value)
	{
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr)
			*ptr = value;
		return *this;
	}

	template<typename T> Vector<T> & Vector<T>::operator += (Vector<T> & v)
	{
		T * sourcePtr = v.Buffer();
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr, ++sourcePtr)
			*ptr += *sourcePtr;
		return *this;
	}

	template<typename T> Vector<T> & Vector<T>::operator -= (Vector<T> & v)
	{
		T * sourcePtr = v.Buffer();
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr, ++sourcePtr)
			*ptr -= *sourcePtr;
		return *this;
	}

	template<typename T> Vector<T> & Vector<T>::operator *= (T value)
	{
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr)
			*ptr *= value;
		return *this;
	}

	template<typename T> Vector<T> & Vector<T>::operator /= (T value)
	{
		for(T * ptr = this->dataPtr; ptr < this->dataEndPtr; ++ptr)
			*ptr /= value;
		return *this;
	}

	template<typename T> Vector<T> operator + (Vector<T> & v1, Vector<T> & v2)
	{
		Vector<T> v = v1;
		v += v2;
		return v;
	}

	template<typename T> Vector<T> operator - (Vector<T> & v1, Vector<T> & v2)
	{
		Vector<T> v = v1;
		v -= v2;
		return v;
	}

	template<typename T> Vector<T> operator + (Vector<T> & vector)
	{
		Vector<T> v = Vector;
		return v;
	}

	template<typename T> Vector<T> operator - (Vector<T> & vector)
	{
		Vector<T> v(Vector.size);
		T * sourcePtr = Vector.dataPtr;
		for(T * ptr = v.dataPtr; ptr < v.dataEndPtr; ++ptr, ++sourcePtr)
			*ptr = -*sourcePtr;
		return v;
	}

	template<typename T> Matrix<T> operator * (Vector<T> & m1, Matrix<T> & m2)
	{
		int c1 = 1;
		int r1 = m1.size;
		int c2 = m2.Col();
		int r2 = m2.Row();
		if (c1 != r2) throw;
		Matrix<T> m(r1, c2);
		int x, y;
		T * m1sizePtr = m1.dataPtr, *ptr1;
		T * m2Ptr = m2.Buffer(), *m2ColPtr, *ptr2;
		T * ptr = m.Buffer();

		for(y = 0; y < r1; ++y, m1sizePtr += c1) //循环m的行
		{
			m2ColPtr = m2Ptr;
			for(x = 0; x < c2; ++x, ++m2ColPtr, ++ptr) //循环m的列
			{
				ptr1 = m1sizePtr;
				ptr2 = m2ColPtr;

				*ptr = *ptr1 * *ptr2;;
			}
		}

		return m;
	}

	template<typename T> Vector<T> operator * (Matrix<T> & m1, Vector<T> & m2)
	{
		int c1 = m1.Col();
		int r1 = m1.Row();
		int r2 = m2.Size();

		if (c1 != r2) throw;
		Vector<T> m(r1);
		int y, z;
		T * m1RowPtr = m1.Buffer(), *ptr1;
		T * m2Ptr = m2.Buffer(), *ptr2;
		T * ptr = m.dataPtr;
		T sum;
		for(y = 0; y < r1; ++y, m1RowPtr += c1, ++ptr) //循环m的行
		{
			ptr1 = m1RowPtr;
			ptr2 = m2Ptr;
			sum = 0;
			for(z = 0; z < c1; ++z, ++ptr1, ++ptr2) //计算m每个单元格的值
			{
				if (*ptr1 * *ptr2 < 0)
				{
					sum = sum;
				}
				sum += *ptr1 * *ptr2;
			}
			*ptr = sum;

		}

		return m;
	}

	template<typename T> T operator * (Vector<T> & v1, Vector<T> & v2)
	{
		if (v1.Size() != v2.Size()) throw;

		T * p1, * p2;
		T sum = 0;
		enumArray2(v1.Buffer(), v2.Buffer(), 0, v1.Size(), p1, p2)
			sum += *p1 * *p2;
		enumArrayEnd;

		return sum;
	}

	template<typename T> T operator * (VectorH<T> & v1, Vector<T> & v2)
	{
		if (v1.Size() != v2.Size()) throw;

		T * p1, * p2;
		T sum = 0;
		enumArray2(v1.Buffer(), v2.Buffer(), 0, v1.Size(), p1, p2)
			sum += *p1 * *p2;
		enumArrayEnd;

		return sum;
	}

	template<typename T> Matrix<T> operator * (Vector<T> & v1, VectorH<T> & v2)
	{
		Matrix<T> m(v1.Size(), v2.Size());

		T * p1, * p2, * p;
		T * p1E = v1.Buffer() + v1.Size();
		T * p2E = v2.Buffer() + v2.Size();

		p = m.Buffer();
		for(p1 = v1.Buffer(); p1 < p1E; ++p1)
			for(p2 = v2.Buffer(); p2 < p2E; ++p2)
				*(p++) = *p1 * *p2;


		return m;
	}

	template<typename T> Vector<T> operator * (T val, Vector<T> & vector)
	{
		Vector<T> m = vector;
		m *= val;

		return m;
	}

	template<typename T> Vector<T> operator * (Vector<T> & vector, T val)
	{
		Vector<T> m = vector;
		m *= val;

		return m;
	}

	template<typename T> Vector<T> operator / (Vector<T> & vector, T val)
	{
		Vector<T> m = vector;
		m /= val;

		return m;
	}

	template<typename T> T & Vector<T>::operator[] (int index)
	{
		return this->dataPtr[index];
	}

	//获得单位向量，该项亮oneIndex处为1，其余为0。size为向量长度，oneIndex为第几个向量值为1（从0开始）
	template<typename T> Vector<T> UnitVector(int size, int oneIndex)
	{
		Vector<T> v(size, (T)0);
		v[oneIndex] = 1;
		return v;
	}

	//获得一个随机向量
	template<typename T> Vector<T> RandomVector(int size, T minValue, T maxValue)
	{
		Vector<T> v(size);
		T * ptr;
		enumArray(v.Buffer(), 0, size, ptr) 
			*ptr = rand(minValue, maxValue);
		enumArrayEnd;

		return v;
	}

	//获得一个随机的单位向量
	template<typename T> Vector<T> RandomUnitVector(int size)
	{
		if (T(1) / 2 == 0)
			return UnitVector<T>(size, rand<int>(0, size));

		Vector<T> v(size);
		T * ptr;
		enumArray(v.Buffer(), 0, size, ptr) 
			*ptr = (T)rand(-1.0, 1.0);
		enumArrayEnd;

		return v / v.Modulo();
	}

	//从向量中截取一部分成为新的向量
	template<typename T>
	Vector<T> Vector<T>::Cut(int start, int len)
	{
		Vector<T> v(len);
		memcpy(v.dataPtr, this->dataPtr + start, len * sizeof(T));
		return v;
	}

	//从向量中截取一部分成为新的向量
	template<typename T>
	Vector<T> Vector<T>::Cut(int start)
	{
		int len = this->size - start;
		Vector<T> v(len);
		memcpy(v.dataPtr, this->dataPtr + start, len * sizeof(T));
		return v;
	}

	//将两个向量合并为一个向量
	template<typename T> Vector<T> MergeVector(Vector<T> & m1, Vector<T> & m2)
	{
		int c1 = m1.Size();
		int c2 = m2.Size();
		int c = c1 + c2;

		Vector<T> m(c, 0);
		
		memcpy(m.Buffer(), m1.Buffer(), m1.MemorySize());
		memcpy((byte*)m.Buffer() + m1.MemorySize(), m2.Buffer(), m2.MemorySize());

		return m;
	}
}