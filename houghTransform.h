#pragma once
#include <cstdlib>
#include <convert.h>
#include <common.h>
#include "memory_.h"

using namespace std;


namespace xmv
{
	struct houghLine;

	template <typename T>
	void houghTransform_CreateMatrix(T * edge, int width, int height, int x1, int x2, int y1, int y2,
		double a1, double a2, double aStep, int stepx, int stepy, int aSize, 
		int dSize, int dMid, double *arr, double *sinArr, double * cosArr, float * directionArr = NULL)
	{
		int aIndex, dIndex;
		double a;

		::memset(arr, 0, aSize * dSize * sizeof(double));
		T* ptrEdgeLine = edge + y1 * width;
		T* ptrEdge;
		float * ptrDLine, * ptrD;
		if (directionArr)
		{
			ptrDLine = directionArr + y1 * width;
		}
		float direction;
		double ins;
		int dx, dy, by, bx, id;
		for (by = y1; by < y2; by += stepy, ptrEdgeLine += width * stepy, directionArr && (ptrDLine += width * stepy))
		{
			ptrEdge = ptrEdgeLine + x1;
			if (directionArr) ptrD = ptrDLine + x1;
			for(bx = x1; bx < x2; bx += stepx, ptrEdge += stepx,  (directionArr && (ptrD += stepx)))
			{
				ins = *ptrEdge;
				if (directionArr)
				{
					direction = *ptrD;
					if (direction > 0)
					{
						direction = direction;
					}
				}

				if (ins <= 0) continue;

				aIndex = 0;
				for(a = a1; aIndex < aSize; a += aStep, ++aIndex)
				{
					if (directionArr)
						if (Abs<T>(sin(a + direction)) > 0.3 && Abs<T>(sin(a - direction)) > 0.3) continue;
					dx = bx - x1; 
					dy = by - y1;

					dIndex = (int)(dx * sin(a) - dy * cos(a) + dMid);
					if (dIndex < 0) continue;
					if (dIndex >= dSize) continue;
					id = aIndex * dSize + dIndex;

					arr[id] += ins;
					if (id > 1) arr[id-1] += ins * 0.8;
					if (id < aSize * dSize - 1) arr[id+1] += ins * 0.8;
					if (id > 2) arr[id-2] += ins * 0.5;
					if (id < aSize * dSize - 2) arr[id+2] += ins * 0.5;
				}
			}
		}
	}

	template <typename T>
	houghLine houghTransform_FindLine(T * edge, int width, int height, int x1, int x2, int y1, int y2,
		double a1, double a2, double aStep, int stepx, int stepy, int aSize, 
		int dSize, int dMid, double *arr, double *sinArr, double * cosArr, int & aMax, int & dMax)
	{
		houghLine l;
		int aIndex, dIndex;

		//findMax
		dMax = 0;
		aMax = 0;
		double maxValue = *arr;
		double *ptr = arr;
		//	CxImage ti(dIndex, aIndex, 8, CXIMAGE_FORMAT_BMP);
		//byte*tibmp = ti.GetData();
		for(aIndex = 0; aIndex < aSize; ++aIndex)
		{
			for(dIndex = 0; dIndex < dSize; ++dIndex, ++ptr)
			{
				/*		int value = *ptr / 10;
				if (value > 255) value = 255;
				if (value < 0) value = 0;
				*tibmp = value;*/
				if (maxValue < *ptr)
				{
					maxValue = *ptr;
					dMax = dIndex;
					aMax = aIndex;
				}
			}
		}
		//	ti.Save("D:\\1.bmp", CXIMAGE_FORMAT_BMP);

		if (fabs(sinArr[aMax]) > 0.7)
		{
			l.x1 = x1 + (dMax - dMid) / sinArr[aMax];
			l.y1 = y1;
			l.x2 = l.x1 + (y2 - y1) * cosArr[aMax] / sinArr[aMax];
			l.y2 = y2;
		}
		else
		{
			l.y1 = y1 - (dMax - dMid) / cosArr[aMax];
			l.x1 = x1;
			l.y2 = l.y1 + (x2 - x1) * sinArr[aMax] / cosArr[aMax];
			l.x2 = x2;
		}
		l.ins = maxValue;

		return l;
	}

	template <typename T>
	houghLine houghTransform(T * edge, int width, int height, int x1, int x2, int y1, int y2,
		double a1, double a2, double aStep, int stepx = 2, int stepy = 2)
	{


		int aSize = (int)((a2 - a1) / aStep);
		int dSize = (x2 - x1) + (y2 - y1) * 2 + 1;
		int dMid = (x2 - x1) + (y2 - y1);

		int aIndex; 
		//	int dIndex;
		//	double d;

		//	int x, y;

		double * sinArr = NEW("XMV.HoughTransform 1: SinArr") double[aSize];
		double * cosArr = NEW("XMV.HoughTransform 1: CosArr") double[aSize];
		double * arr =    NEW("XMV.HoughTransform 1: arr")	  double[aSize * dSize];

		aIndex = 0;
		double a;
		for(a = a1; aIndex < aSize; a += aStep, ++aIndex)
		{
			sinArr[aIndex] = sin(a);
			cosArr[aIndex] = cos(a);
		}

		int dMax, aMax;

		houghTransform_CreateMatrix(edge, width, height, x1, x2, y1, y2, a1, a2, aStep, stepx, stepy, 
			aSize, dSize, dMid, arr, sinArr, cosArr);

		houghLine l = houghTransform_FindLine(edge, width, height, x1, x2, y1, y2, a1, a2, aStep, stepx, stepy, 
			aSize, dSize, dMid, arr, sinArr, cosArr, aMax, dMax);

		DELETEARR(arr);
		DELETEARR(sinArr);
		DELETEARR(cosArr);

		return l;
	}

	template <typename T>
	int houghTransform(T * edge, houghLine * lines, int lineCount, int width, int height, int x1, int x2, int y1, int y2,
		double a1, double a2, double aStep, double deepth = 0.2, int stepx = 2, int stepy = 2, bool absDepth = false, 
		float * directionArr = NULL)
	{
		int aSize = (int)((a2 - a1) / aStep);
		int dSize = (x2 - x1) + (y2 - y1) * 2 + 1;
		int dMid = (x2 - x1) + (y2 - y1);

		int aIndex; 
		//	int dIndex;
		//	double d;

		//	int y;

		double * sinArr = NEW("XMV.HoughTransform 2: SinArr")  double[aSize];
		double * cosArr = NEW("XMV.HoughTransform 2: CosArr")  double[aSize];
		double * arr    = NEW("XMV.HoughTransform 2: arr")	   double[aSize * dSize];

		aIndex = 0;
		double a;
		for(a = a1; aIndex < aSize; a += aStep, ++aIndex)
		{
			sinArr[aIndex] = sin(a);
			cosArr[aIndex] = cos(a);
		}

		int dMax, aMax;

		houghTransform_CreateMatrix(edge, width, height, x1, x2, y1, y2, a1, a2, aStep, stepx, stepy, 
			aSize, dSize, dMid, arr, sinArr, cosArr, directionArr);

		double val = Max(arr, dSize, aSize, 0, dSize, 0, aSize, 2, 2) * deepth;
		if (absDepth) val = deepth ;

		int c = 0;
		for(int i = 0; i < lineCount; ++i, ++c)
		{
			houghLine l = houghTransform_FindLine(edge, width, height, x1, x2, y1, y2, a1, a2, aStep, stepx, stepy, 
				aSize, dSize, dMid, arr, sinArr, cosArr, aMax, dMax);

			if (arr[aMax * dSize + dMax] < val) break;

			*lines++ = l;
			::ClearLightArea(arr, dSize, aSize, 0, dSize, 0, aSize, dMax, aMax, val);
		}

		DELETEARR(arr);
		DELETEARR(sinArr);
		DELETEARR(cosArr);

		return c;
	}

}