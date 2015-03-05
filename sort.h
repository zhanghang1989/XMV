
#pragma once

#include <stdlib.h>



namespace xmv
{
	template<typename T>
	void BubbleSort(T * arr, int size)
	{
		T * end = arr + size;
		T * p1 = arr;
		T * p2 = p1 + 1;
		T temp;

		bool exchanged = false;
		do
		{
			p1 = arr;
			p2 = p1 + 1;
			exchanged = false;
			for(; p2 < end; ++p1, ++p2)
			{
				if (*p1 > *p2)
				{
					temp = *p1;
					*p1 = *p2;
					*p2 = temp;
					exchanged = true;
				}
			}
		}
		while(exchanged);
	}

	template<typename T>
	void BubbleSort(T * arr, int size, int compare(T & x, T & y))
	{
		T * end = arr + size;
		T * p1 = arr;
		T * p2 = p1 + 1;
		T temp;

		bool exchanged = false;
		do
		{
			p1 = arr;
			p2 = p1 + 1;
			exchanged = false;
			for(; p2 < end; ++p1, ++p2)
			{
				if (compare(*p1, *p2) > 0)
				{
					temp = *p1;
					*p1 = *p2;
					*p2 = temp;
					exchanged = true;
				}
			}
		}
		while(exchanged);
	}

	template <typename T>
	void Quicksort(T * start, T * end)		//快速排序函数
	{
		if (start >= end - 1) return;				//数组长度小于等于1，无需排序

		T key = *start;						//获取主元
		T temp;								//交换数据时需要的临时变量
		T * left = start + 1;					//left指针指向数组下标1的元素
		T * right = end - 1;						//right指针指向数组的结尾

		for(; ; )								//循环交换
		{
			for( ; left < end; left++)			//left向右寻找大于等于key的数据
				if (*left >= key) break;			//找到大于等于key的数据，跳出循环
			for( ; right > start; right--)		//right向左找到小于key的数据
				if (*right < key) break;			//找到小于key的数据，跳出循环

			if (left >= right) 					//如果left移动到了right的右侧
				break;							//终止循环，本轮交换结束

			temp = *left;						//对两个数字进行交换
			*left = *right;
			*right = temp;

			left++;				//交换后，left向后移动一个，准备下一次查找
			right--;				//交换后，right向前移动一个，准备下一次查找
		}

		temp = * start;							//将主元放入两组数据的中间
		*start = *right;
		*right = temp;

		Quicksort(start, right);				//对主元左侧分组进行下一轮排序
		Quicksort(right + 1, end);				//对主元右侧分组进行下一轮排序
	}

	template <typename T>
	void Quicksort(T * arr, int size)	//函数重载，可用习惯的方法调用排序函数
	{
		Quicksort(arr, arr + size);
	}

	template <typename T>
	void Quicksort(T * start, T * end, int compare(T & x, T & y))		//快速排序函数
	{
		if (start >= end - 1) return;				//数组长度小于等于1，无需排序

		T key = *start;						//获取主元
		T temp;								//交换数据时需要的临时变量
		T * left = start + 1;					//left指针指向数组下标1的元素
		T * right = end - 1;						//right指针指向数组的结尾

		for(; ; )								//循环交换
		{
			for( ; left < end; left++)			//left向右寻找大于等于key的数据
				if (compare(*left, key) > 0) break;			//找到大于等于key的数据，跳出循环
			for( ; right > start; right--)		//right向左找到小于key的数据
				if (compare(*right, key) < 0) break;			//找到小于key的数据，跳出循环

			if (left >= right) 					//如果left移动到了right的右侧
				break;							//终止循环，本轮交换结束

			temp = *left;						//对两个数字进行交换
			*left = *right;
			*right = temp;

			left++;				//交换后，left向后移动一个，准备下一次查找
			right--;				//交换后，right向前移动一个，准备下一次查找
		}

		temp = * start;							//将主元放入两组数据的中间
		*start = *right;
		*right = temp;

		Quicksort(start, right, compare);				//对主元左侧分组进行下一轮排序
		Quicksort(right + 1, end, compare);				//对主元右侧分组进行下一轮排序
	}

	template <typename T>
	void Quicksort(T * arr, int size, int compare(T & x, T & y))
	{
		Quicksort(arr, arr + size, compare);
	}

	template <typename T>
	int Find(T * arr, int size,  T & value)
	{
		T * p;
		enumArray(T * arr, 0, size, p)
			if (*p == value) return p - arr;
		enumArrayEnd;

		return -1;
	}

	template <typename T>
	int DichotomyFind(T * arr, int start, int end, T & value)
	{
		if (end <= start) return -1;
		if (arr[start] == value) return start;
		if (arr[end - 1] == value) return end;
		if (end - start <= 2) return -1;
		int mid = (start + end - 1) >> 1;
		T & v = arr[mid];
		if (value == v) return mid;
		else if (value < v) return DichotomyFind(arr, start, mid, value);
		else return DichotomyFind(arr, mid + 1, end, value);
	}

	template <typename T>
	int DichotomyFind(T * arr, int size, T & value)
	{
		return DichotomyFind(arr, 0, size, value);
	}
}