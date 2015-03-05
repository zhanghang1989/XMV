#pragma once

#include <cstdlib>
#include "pixcel.h"
#include "memory.h"
#include <ctime>

using namespace std;

typedef unsigned char byte;

#define CLR

namespace xmv
{


#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#define checkBound(bx1,bx2,by1,by2,x1,x2,y1,y2) \
	if ((x1) < (bx1)) (x1) = (bx1); \
	if ((x1) > (bx2)) (x1) = (bx2); \
	if ((x2) < (bx1)) (x2) = (bx1); \
	if ((x2) > (bx2)) (x2) = (bx2); \
	if ((y1) < (by1)) (y1) = (by1); \
	if ((y1) > (by2)) (y1) = (by2); \
	if ((y2) < (by1)) (y2) = (by1); \
	if ((y2) > (by2)) (y2) = (by2); 
#define drawPoint(bmp,imgWidth,imgHeight,x,y,value) (bmp[(y) * (imgWidth) + (x)] = (value));
#define between(v,v1,v2) ((v) >= (v1) && (v) <= (v2))
#define EPS (1E-10)

	template<typename T> inline void ClearImage(T * img, int size)
	{
		::memset(img, 0, size  * sizeof(T));
	}

	template<typename T> inline void ClearImage(T * img, int width, int height)
	{
		::memset(img, 0, width * height * sizeof(T));
	}

	template<typename T> void ClearImage(T * img, int width, int height, int x1, int x2, int y1, int y2)
	{
		checkBound(0, width, 0, height, x1, x2, y1, y2);
		int y;
		T * ptrLine = img + y1 * width + x1;
		int len = (x2 - x1) * sizeof(T);
		for(y = y1; y < y2; ++y, ptrLine += width)
		{
			::memset(ptrLine, 0, len);
		}
	}



	template<typename T, typename T2> void CopyImage(T * targetImg, int size, T2 * source)
	{
		T * ptr = targetImg;
		T * ptrEnd = targetImg + size;
		T2 * t = source;

		for(; ptr < ptrEnd; ++ptr, ++t)
		{
			*ptr = (T)*t;
		}
	}

	struct houghLine
	{
		double x1;double x2;double y1;double y2;
		houghLine * left;
		houghLine * right;
		double ins;
	};


	template<typename T> void DrawLine(T *bmp, int imgWidth, int imgHeight, houghLine & l, T value)
	{
		DrawLine(bmp, imgWidth, imgHeight, l.x1, l.y1, l.x1, l.y2, value);
	}

	template<typename T, typename T2> void CopyAndAddImage(T * target, int targetWidth, int targetHeight,
		int targetX, int targetY, T2 * source, int sourceWidth, int sourceHeight)
	{
		int y, targetYEnd = targetY + sourceHeight;
		int x, targetXEnd = targetX + sourceWidth;
		T  * targetLinePtr = target + targetY * targetWidth + targetX, *t;
		T2 * s = source;

		for(y = targetY; y < targetYEnd; ++y, targetLinePtr += targetWidth)
		{
			if (y < 0) { s += sourceWidth; continue; }
			if (y >= targetHeight) { s += sourceWidth; continue; }

			t = targetLinePtr;
			for(x = targetX; x < targetXEnd; ++x, ++t, ++s)
			{
				if (x < 0) continue;
				if (x >= targetWidth) continue;
				*t += *s;
			}
		}
	}

	template<typename T, typename T2> void CopyImage(T * target, int targetWidth, int targetHeight,
		int targetX, int targetY, T2 * source, int sourceWidth, int sourceHeight)
	{
		int y, targetYEnd = targetY + sourceHeight;
		int x, targetXEnd = targetX + sourceWidth;
		T  * targetLinePtr = target + targetY * targetWidth + targetX, *t;
		T2 * s = source;

		for(y = targetY; y < targetYEnd; ++y, targetLinePtr += targetWidth)
		{
			if (y < 0) { s += sourceWidth; continue; }
			if (y >= targetHeight) { s += sourceWidth; continue; }

			t = targetLinePtr;
			for(x = targetX; x < targetXEnd; ++x, ++t, ++s)
			{
				if (x < 0) continue;
				if (x >= targetWidth) continue;
				*t = *s;
			}
		}
	}


	template<typename T, typename T2> void CopyImage(T * target, int width, int height, int x1, int x2, int y1, int y2, T2 * source)
	{
		CopyImage(target, width, height, x1, y1, source, x2 - x1, y2 - y1);
	}


	template<typename T, typename T2> void SetImageValue(T * img, int width, int height, int x1, int x2, int y1, int y2, T2 value)
	{
		int y;
		int x;
		T * ptrLine = img + y1 * width + x1;
		T * ptr;

		for(y = y1; y < y2; ++y, ptrLine += width)
		{
			ptr = ptrLine;
			for(x = x1; x < x2; ++x, ++ptr)
				*ptr = (T)value;
		}
	}

	template<typename T> int GetImageValueCount(T * img, int width, int height, int x1, int x2, int y1, int y2)
	{
		int y;
		int x;
		int c = 0;
		T * ptrLine = img + y1 * width + x1;
		T * ptr;

		for(y = y1; y < y2; ++y, ptrLine += width)
		{
			ptr = ptrLine;
			for(x = x1; x < x2; ++x, ++ptr)
				if (*ptr) c++;
		}
		return c;
	}

	template<typename T> void DrawHorizontalLine(T * img, int width, int height, int x1, int x2, int y, T value)
	{
		T * ptr = img + y * width + x1;
		T * end = img + y * width + x2;

		for(; ptr < end; ++ptr)
			*ptr = value;
	}

	template<typename T> void DrawVerticalLine(T * img, int width, int height, int x, int y1, int y2, T value)
	{
		T * ptr = img + y1 * width + x;
		T * end = img + y2 * width + x;

		for(; ptr < end; ptr += width)
			*ptr = value;
	}

	template<typename T> void DrawCross(T *bmp, int w, int h, int x0, int y0, int size, T value, byte alpha = 255)
	{
		DrawLine<T>(bmp, w, h, x0 - size, y0, x0 + size + 1, y0, value, alpha);
		DrawLine<T>(bmp, w, h, x0, y0 - size, x0, y0 + size + 1, value, alpha);
	}

	template<typename T> void DrawLine(T *bmp, int imgWidth, int imgHeight, int x1, int y1, int x2, int y2, T value, byte alpha = 255)
	{
		double k;

		int x, y;
		if (Abs(x1 - x2) >= Abs(y1 - y2))
		{
			//x based
			k = (double)(y2 - y1) / (x2 - x1);
			if (x1 <= x2)
				for(x = (int)x1; x < x2; x++)
				{
					y = (int)((x - x1) * k + y1 + 0.5F);
					if (x < 0) continue;
					if (x >= imgWidth) continue;
					if (y < 0) continue;
					if (y >= imgHeight) continue;
					bmp[y * imgWidth + x] = value * alpha / 255 + bmp[y * imgWidth + x] * (255 - alpha) / 255;
				}
			else
				for(x = (int)x1; x >= x2; x--)
				{
					y = (int)((x - x1) * k + y1 + 0.5F);
					if (x < 0) continue;
					if (x >= imgWidth) continue;
					if (y < 0) continue;
					if (y >= imgHeight) continue;
					bmp[y * imgWidth + x] = value * alpha / 255 + bmp[y * imgWidth + x] * (255 - alpha) / 255;
				}
		}
		else
		{
			//y based 
			k = (double)(x2 - x1) / (y2 - y1);
			if (y1 <= y2)
				for(y = (int)y1; y < y2; y++)
				{
					x = (int)((y - y1) * k + x1);
					if (x < 0) continue;
					if (x >= imgWidth) continue;
					if (y < 0) continue;
					if (y >= imgHeight) continue;
					bmp[y * imgWidth + x] = value * alpha / 255 + bmp[y * imgWidth + x] * (255 - alpha) / 255;
				}
			else
				for(y = (int)y1; y >= y2; y--)
				{
					x = (int)((y - y1) * k + x1);

					if (x < 0) continue;
					if (x >= imgWidth) continue;
					if (y < 0) continue;
					if (y >= imgHeight) continue;
					bmp[y * imgWidth + x] = value * alpha / 255 + bmp[y * imgWidth + x] * (255 - alpha) / 255;
				}
		}
	}


	template<typename T> void DrawRectangle(T * img, int width, int height, int x1, int x2, int y1, int y2, T value)
	{
		DrawHorizontalLine(img, width, height, x1, x2 + 1, y1, value);
		DrawHorizontalLine(img, width, height, x1, x2 + 1, y2, value);
		DrawVerticalLine(img, width, height, x1, y1, y2 + 1, value);
		DrawVerticalLine(img, width, height, x2, y1, y2 + 1, value);
	}

	template<typename T> void FillRectangle(T * img, int width, int height, int x1, int x2, int y1, int y2, T value)
	{
		T * p;
		enumImage(img, width, height, x1, x2, y1, y2, p)
			*p = value;
		enumImageEnd;
	}


	extern bool ____randomizeOk;
	template<typename T> inline T rand(T max)
	{
		if (! ____randomizeOk) { std::srand(time(0)); std::rand(); std::rand(); ____randomizeOk = true; }
		return (T)(std::rand() * max / RAND_MAX);
	}

	template<typename T> inline T rand(T min, T max)
	{
		if (! ____randomizeOk) { std::srand(time(0)); std::rand(); std::rand(); ____randomizeOk = true; }
		return (T)(std::rand() * (max - min) / RAND_MAX) + min;
	}
}

