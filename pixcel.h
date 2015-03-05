#pragma once;

#include "common.h"

using namespace std;
namespace xmv
{
#define enumImageEnd }}}

#define enumImage(bmp, width, height, x1, x2, y1, y2, pointer)  \
	{ \
	int x_____, y_____; \
	int x1_____ = (x1); \
	int y1_____ = (y1); \
	int x2_____ = (x2); \
	int y2_____ = (y2); \
	int stepy_____ = (width)-(x2_____ - x1_____); \
	int offset_____ = y1_____ * width + x1_____; \
	pointer = (bmp)+offset_____; \
	for (y_____ = y1_____; y_____ < y2_____; ++y_____, pointer += stepy_____) \
	{ \
	for (x_____ = x1_____; x_____ < x2_____; ++x_____, ++pointer) \
	{

#define enumImageStep(bmp, width, height, x1, x2, y1, y2, xStep, yStep, pointer)  \
	{ \
	int x_____, y_____; \
	int x1_____ = (x1); \
	int y1_____ = (y1); \
	int x2_____ = (x2); \
	int y2_____ = (y2); \
	int xs_____ = (xStep); \
	int ys_____ = (yStep); \
	int stepy_____ = ys_____ * (width)-(x2_____ - x1_____) / xs_____ * xs_____; \
	int offset_____ = y1_____ * (width)+x1_____; \
	pointer = (bmp)+offset_____; \
	for (y_____ = y1_____; y_____ < y2_____; y_____ += ys_____, pointer += stepy_____) \
	{ \
	for (x_____ = x1_____; x_____ < x2_____; x_____ += xs_____, pointer += xs_____) \
	{

#define enumImage2(bmp1, bmp2, width, height, x1, x2, y1, y2, pointer1, pointer2)  \
	{ \
	int x_____, y_____; \
	int x1_____ = (x1); \
	int y1_____ = (y1); \
	int x2_____ = (x2); \
	int y2_____ = (y2); \
	int stepy_____ = (width)-(x2_____ - x1_____); \
	int offset_____ = y1_____ * (width)+x1_____; \
	pointer1 = (bmp1)+offset_____; \
	pointer2 = (bmp2)+offset_____; \
	for (y_____ = y1_____; y_____ < y2_____; ++y_____, pointer1 += stepy_____, pointer2 += stepy_____) \
	{ \
	for (x_____ = x1_____; x_____ < x2_____; ++x_____, ++pointer1, ++pointer2) \
	{

#define enumImage3(bmp1, bmp2, bmp3, width, height, x1, x2, y1, y2, pointer1, pointer2, pointer3)  \
	{ \
	int x_____, y_____; \
	int x1_____ = (x1); \
	int y1_____ = (y1); \
	int x2_____ = (x2); \
	int y2_____ = (y2); \
	int stepy_____ = (width)-(x2_____ - x1_____); \
	int offset_____ = y1_____ * (width)+x1_____; \
	pointer1 = (bmp1)+offset_____; \
	pointer2 = (bmp2)+offset_____; \
	pointer3 = (bmp3)+offset_____; \
	for (y_____ = y1_____; y_____ < y2_____; ++y_____, \
	pointer1 += stepy_____, \
	pointer2 += stepy_____, \
	pointer3 += stepy_____) \
	{ \
	for (x_____ = x1_____; x_____ < x2_____; ++x_____, \
	++pointer1, \
	++pointer2, \
	++pointer3) \
	{


#define enumImage4(bmp1, bmp2, bmp3, bmp4, width, height, x1, x2, y1, y2, pointer1, pointer2, pointer3, pointer4)  \
	{ \
	int x_____, y_____; \
	int x1_____ = (x1); \
	int y1_____ = (y1); \
	int x2_____ = (x2); \
	int y2_____ = (y2); \
	int stepy_____ = (width)-(x2_____ - x1_____); \
	int offset_____ = y1_____ * (width)+x1_____; \
	pointer1 = (bmp1)+offset_____; \
	pointer2 = (bmp2)+offset_____; \
	pointer3 = (bmp3)+offset_____; \
	pointer4 = (bmp4)+offset_____; \
	for (y_____ = y1_____; y_____ < y2_____; ++y_____, \
	pointer1 += stepy_____, \
	pointer2 += stepy_____, \
	pointer3 += stepy_____, \
	pointer4 += stepy_____) \
	{ \
	for (x_____ = x1_____; x_____ < x2_____; ++x_____, \
	++pointer1, \
	++pointer2, \
	++pointer3, \
	++pointer4) \
	{


#define enumImage5(bmp1, bmp2, bmp3, bmp4, bmp5, width, height, x1, x2, y1, y2, pointer1, pointer2, pointer3, pointer4, pointer5)  \
	{ \
	int x_____, y_____; \
	int x1_____ = (x1); \
	int y1_____ = (y1); \
	int x2_____ = (x2); \
	int y2_____ = (y2); \
	int stepy_____ = (width)-(x2_____ - x1_____); \
	int offset_____ = y1_____ * (width)+x1_____; \
	pointer1 = (bmp1)+offset_____; \
	pointer2 = (bmp2)+offset_____; \
	pointer3 = (bmp3)+offset_____; \
	pointer4 = (bmp4)+offset_____; \
	pointer5 = (bmp5)+offset_____; \
	for (y_____ = y1_____; y_____ < y2_____; ++y_____, \
	pointer1 += stepy_____, \
	pointer2 += stepy_____, \
	pointer3 += stepy_____, \
	pointer4 += stepy_____, \
	pointer5 += stepy_____) \
	{ \
	for (x_____ = x1_____; x_____ < x2_____; ++x_____, \
	++pointer1, \
	++pointer2, \
	++pointer3, \
	++pointer4, \
	++pointer5) \
	{

#define enumImage6(bmp1, bmp2, bmp3, bmp4, bmp5, bmp6, width, height, x1, x2, y1, y2, pointer1, pointer2, pointer3, pointer4, pointer5, pointer6)  \
	{ \
	int x_____, y_____; \
	int x1_____ = (x1); \
	int y1_____ = (y1); \
	int x2_____ = (x2); \
	int y2_____ = (y2); \
	int stepy_____ = (width)-(x2_____ - x1_____); \
	int offset_____ = y1_____ * (width)+x1_____; \
	pointer1 = (bmp1)+offset_____; \
	pointer2 = (bmp2)+offset_____; \
	pointer3 = (bmp3)+offset_____; \
	pointer4 = (bmp4)+offset_____; \
	pointer5 = (bmp5)+offset_____; \
	pointer6 = (bmp6)+offset_____; \
	for (y_____ = y1_____; y_____ < y2_____; ++y_____, \
	pointer1 += stepy_____, \
	pointer2 += stepy_____, \
	pointer3 += stepy_____, \
	pointer4 += stepy_____, \
	pointer5 += stepy_____) \
	pointer6 += stepy_____) \
	{ \
	for (x_____ = x1_____; x_____ < x2_____; ++x_____, \
	++pointer1, \
	++pointer2, \
	++pointer3, \
	++pointer4, \
	++pointer5, \
	++pointer6, \
	) \
	{


#define enumArray(arr, x1, x2, pointer) \
	{ \
	auto end_____ = (arr)+(x2); \
	for ((pointer) = (arr)+(x1); (pointer) < (end_____); ++(pointer)) \
	{

#define enumArray2(arr1, arr2, x1, x2, pointer1, pointer2) \
	{ \
	auto end_____ = (arr1)+(x2); \
	(pointer1) = (arr1)+(x1); \
	(pointer2) = (arr2)+(x1); \
	for (; (pointer1) < (end_____); ++pointer1, ++pointer2) \
	{

#define enumArray3(arr1, arr2, arr3, x1, x2, pointer1, pointer2, pointer3) \
	{ \
	void * end_____ = (arr1)+(x2); \
	(pointer1) = (arr1)+(x1); \
	(pointer2) = (arr2)+(x1); \
	(pointer3) = (arr3)+(x1); \
	for (; (pointer1) < (end_____); ++(pointer1), ++(pointer2), ++(pointer3)) \
	{

#define enumArray4(arr1, arr2, arr3, arr4, x1, x2, pointer1, pointer2, pointer3, pointer4) \
	{ \
	void * end_____ = (arr1)+(x2); \
	(pointer1) = (arr1)+(x1); \
	(pointer2) = (arr2)+(x1); \
	(pointer3) = (arr3)+(x1); \
	(pointer4) = (arr4)+(x1); \
	for (; (pointer1) < (end_____); ++(pointer1), ++(pointer2), ++(pointer3), ++(pointer4)) \
	{

#define enumArray5(arr1, arr2, arr3, arr4, arr5, x1, x2, pointer1, pointer2, pointer3, pointer4, pointer5) \
	{ \
	void * end_____ = (arr1)+(x2); \
	(pointer1) = (arr1)+(x1); \
	(pointer2) = (arr2)+(x1); \
	(pointer3) = (arr3)+(x1); \
	(pointer4) = (arr4)+(x1); \
	(pointer5) = (arr5)+(x1); \
	for (; (pointer1) < (end_____); ++(pointer1), ++(pointer2), ++(pointer3), ++(pointer4), ++(pointer5)) \
	{

#define enumArray6(arr1, arr2, arr3, arr4, arr5, arr6, x1, x2, pointer1, pointer2, pointer3, pointer4, pointer5, pointer6) \
	{ \
	void * end_____ = (arr1)+(x2); \
	(pointer1) = (arr1)+(x1); \
	(pointer2) = (arr2)+(x1); \
	(pointer3) = (arr3)+(x1); \
	(pointer4) = (arr4)+(x1); \
	(pointer5) = (arr5)+(x1); \
	(pointer6) = (arr6)+(x1); \
	for (; (pointer1) < (end_____); ++(pointer1), ++(pointer2), ++(pointer3), ++(pointer4), ++(pointer5), ++(pointer6)) \
	{


#define enumArray7(arr1, arr2, arr3, arr4, arr5, arr6, arr7, x1, x2, pointer1, pointer2, pointer3, pointer4, pointer5, pointer6, pointer7) \
	{ \
	void * end_____ = (arr1)+(x2); \
	(pointer1) = (arr1)+(x1); \
	(pointer2) = (arr2)+(x1); \
	(pointer3) = (arr3)+(x1); \
	(pointer4) = (arr4)+(x1); \
	(pointer5) = (arr5)+(x1); \
	(pointer6) = (arr6)+(x1); \
	(pointer7) = (arr7)+(x1); \
	for (; (pointer1) < (end_____); ++(pointer1), ++(pointer2), ++(pointer3), ++(pointer4), ++(pointer5), ++(pointer6), ++(pointer7)) \
	{
#define enumArrayEnd }}


#define rotate_(xFrom, yFrom, alpha, xTo, yTo) \
	(xTo) = (xFrom)* cos(alpha) - (yFrom)* sin(alpha); \
	(yTo) = (xFrom)* sin(alpha) + (yFrom)* cos(alpha);

	//从img图像中选取(x,y)点的值
	template<typename T> inline T Sample(T * img, int width, int height, int x, int y)
	{
		return img[y * width + x];
	}

	//从img图像中选取(x,y)点的值
	template<typename T> inline T Sample(T * img, int width, int height, double dx, double dy)
	{
		if (dx < 0) dx = 0;
		if (dy < 0) dy = 9;
		if (dx >= width - 1) dx = width - 1;
		if (dy >= height - 1) dy = height - 1;


		int Y = dy;
		int X = dx;


		double ky = dy - Y;
		double kx = dx - X;
		if (ky == 0) Y--;
		if (kx == 0) X--;

		int pos = Y * width + X;

		if (pos < 0) return 0;
		if (pos >= width * height) return 0;
		return img[pos] * (1 - kx) * (1 - ky) + img[pos + 1] * kx * (1 - ky)
			+ img[pos + width] * (1 - kx) * ky + img[pos + width + 1] * kx * ky;
	}


	template<typename T> void ClearLightArea(T * img, int width, int height, int x1, int x2, int y1, int y2, int x0, int y0, T value)
	{
		int delXCount = 0;
		T * linePtr;
		T * ptr;
		linePtr = img + width * y0;
		for (int y = y0; y < y2; ++y, linePtr += width)
		{
			delXCount = 0;

			ptr = linePtr + x0;
			for (int x = x0; x >= x1; --x, --ptr)
			{
				if (*ptr > value)
				{
					*ptr = 0;
					delXCount++;
				}
				else break;
			}
			ptr = linePtr + x0 + 1;
			for (int x = x0 + 1; x < x2; ++x, ++ptr)
			{
				if (*ptr >= value)
				{
					*ptr = 0;
					delXCount++;
				}
				else break;
			}
			if (delXCount <= 0) break;
		}

		linePtr = img + width * (y0 - 1);
		for (int y = y0 - 1; y >= y1; --y, linePtr -= width)
		{
			delXCount = 0;

			ptr = linePtr + x0;
			for (int x = x0; x >= x1; --x, --ptr)
			{
				if (*ptr > value)
				{
					*ptr = 0;
					delXCount++;
				}
				else break;
			}
			ptr = linePtr + x0 + 1;
			for (int x = x0 + 1; x < x2; ++x, ++ptr)
			{
				if (*ptr >= value)
				{
					*ptr = 0;
					delXCount++;
				}
				else break;
			}
			if (delXCount <= 0) break;
		}
	}

	template<typename TImg, typename TGH>
	void GetHistogram(TImg *bmp, int width, int height, TImg bmpMaxValue, TGH * histogram, int histogramSize, int x1, int x2, int y1, int y2)
	{
		::ClearImage(histogram, histogramSize);

		TImg * p;
		enumImage(bmp, width, height, x1, x2, y1, y2, p)
			++histogram[*p * (histogramSize - 1) / bmpMaxValue];
		enumImageEnd
	}

	template<typename TImg, typename TI> void HorizontalIntegral(TImg * bmp, int width, int height, int x1, int x2, int y1, int y2, TI * integral, bool divSize = true)
	{
		int hSpan = y2 - y1;
		int wSpan = x2 - x1;

		::ClearImage(integral, hSpan);

		TImg * line = bmp + x1;
		TImg * p;
		TI * ip = integral;
		int x, y;
		for (y = y1; y < y2; ++y, line += width, ++ip)
		{
			p = line;
			for (x = x1; x < x2; ++x, ++p)
			{
				*ip += *p;
			}
		}
		ip = integral;
		TI* ipEnd = ip + hSpan;

		if (divSize)
		for (; ip < ipEnd; ++ip)
			*ip /= wSpan;
	}


	template<typename TImg, typename TI> void VerticalIntegral(TImg * bmp, int width, int height, int x1, int x2, int y1, int y2, TI * integral, bool divSize = true)
	{
		int hSpan = y2 - y1;
		int wSpan = x2 - x1;

		::ClearImage(integral, wSpan);

		TImg * line = bmp + x1;
		TImg * p;
		TI * ip = integral;
		int x, y;
		for (y = y1; y < y2; ++y, line += width)
		{
			p = line;
			ip = integral;
			for (x = x1; x < x2; ++x, ++p, ++ip)
			{
				*ip += *p;
			}
		}
		ip = integral;
		TI* ipEnd = ip + wSpan;

		if (divSize)
		for (; ip < ipEnd; ++ip)
			*ip /= hSpan;
	}

	template <typename T>
	T Max(T *image, int width, int height, int x1, int x2, int y1, int y2, int stepX, int stepY)
	{
		T maxValue = 0;
		T *np, *npy;      // 图像指针

		int i, j;          // various counters
		int lineStep = stepY * width;

		npy = image + width * y1 + x1;
		for (i = y1; i < y2; i += stepY, npy += lineStep)
		{
			np = npy;
			for (j = x1; j < x2; j += stepX, np += stepX)
			{
				if (*np > maxValue) maxValue = *np;
			}
		}
		if (maxValue <= 0) return 0;
		return maxValue;
	}

	template <typename T>
	T Otsu(T *image, int width, int height, int x1, int x2, int y1, int y2, int stepX = 1, int stepY = 1)
	{
		T maxValue = 0;
		T *np, *npy;      // 图像指针
		int thresholdValue = 1; // 阈值
		int ihist[256];             // 图像直方图，256个点

		int i, j, k;          // various counters
		double n, n1, n2, gmin, gmax;
		double m1, m2, sum, csum, fmax, sb;
		int lineStep = stepY * width;

		// 对直方图置零...
		memset(ihist, 0, sizeof(ihist));

		npy = image + width * y1 + x1;
		for (i = y1; i < y2; i += stepY, npy += lineStep)
		{
			np = npy;
			for (j = x1; j < x2; j += stepX, np += stepX)
			{
				if (*np > maxValue) maxValue = *np;
			}
		}
		if (maxValue <= 0) return 0;

		gmin = 255; gmax = 0;
		// 生成直方图
		npy = image + width * y1 + x1;
		int value;
		for (i = y1; i < y2; i += stepY, npy += lineStep)
		{
			np = npy;
			for (j = x1; j < x2; j += stepX, np += stepX)
			{
				value = (int)(*np * 255 / maxValue);
				if (value < 0) value = 0;
				if (value > 255) value = 255;
				if (value > 0)
					ihist[value]++;
				if (value > gmax) gmax = *np;
				if (value < gmin) gmin = *np;
				//np++; /* next pixel */
			}
		}

		// set up everything
		sum = csum = 0.0;
		n = 0;

		for (k = 0; k <= 255; k++) {
			sum += (double)k * (double)ihist[k];    /* x*f(x) 质量矩*/
			n += ihist[k];                            /*  f(x)    质量    */
		}

		if (!n)
		{
			// if n has no value, there is problems...
			fprintf(stderr, "NOT NORMAL thresholdValue = 160\n");
			return (60);
		}

		// do the otsu global thresholding method
		fmax = -1.0;
		n1 = 0;
		for (k = 0; k < 255; ++k)
		{
			n1 += ihist[k];
			if (!n1) { continue; }
			n2 = n - n1;
			if (n2 == 0) { break; }
			csum += (double)k *ihist[k];
			m1 = csum / n1;
			m2 = (sum - csum) / n2;
			sb = (double)n1 *(double)n2 *(m1 - m2) * (m1 - m2);
			/* bbg: note: can be optimized. */
			if (sb > fmax)
			{
				fmax = sb;
				thresholdValue = k;
			}
		}

		// at this point we have our thresholding value

		// debug code to display thresholding values


		return (thresholdValue * maxValue / 255);
	}

	//平滑
	template<typename T>
	void Smooth(T * bmp, int width, int height, int x1, int x2, int y1, int y2, double centerPower = 0.5)
	{
		T * b;
		double p = (1 - centerPower) / 4;
		enumImage(bmp, width, height, x1 + 1, x2 - 1, y1 + 1, y2 - 1, b)
			*b = (T)(*b * centerPower + (b[1] + b[-1] + b[width] + b[-width]) * p);
		enumImageEnd
	}

	//中值滤波
	template<typename T>
	void SmoothMidCorss(T * bmp, int width, int height, int x1, int x2, int y1, int y2)
	{
		T bs[5];
		T * b;
		enumImage(bmp, width, height, x1 + 1, x2 - 1, y1 + 1, y2 - 1, b)
			bs[0] = *b;
		bs[1] = b[1];
		bs[2] = b[-1];
		bs[3] = b[width];
		bs[4] = b[-width];
		BubbleSort(bs, 5);
		*b = bs[2];
		enumImageEnd
	}

	//缩小二值图像，0为背景，1为前景
	template<typename T>
	void Shrink(T * bmp, T * output, int width, int height)
	{
		T * p1, *p2;
		int size = width * height;
		memcpy(output, bmp, sizeof(T)* size);
		enumImage2(bmp, output, width, height, 0, width, 0, height, p1, p2)
		{
			if (!*p1)
			{
				if (x_____ > 0) p2[-1] = 0;
				if (x_____ < width - 1) p2[1] = 0;
				if (y_____ > 0) p2[-width] = 0;
				if (y_____ < height - 1) p2[width] = 0;
			}
		}
		enumImageEnd;
	}

	//扩展二值图像，0为背景，1为前景
	template<typename T>
	void Expand(T * bmp, T * output, int width, int height)
	{
		T * p1, *p2;
		int size = width * height;
		memcpy(output, bmp, sizeof(T)* size);
		enumImage2(bmp, output, width, height, 0, width, 0, height, p1, p2)
		{
			if (*p1 > 0.0)
			{
				if (x_____ > 0) p2[-1] = *p1;
				if (x_____ < width - 1) p2[1] = *p1;
				if (y_____ > 0) p2[-width] = *p1;
				if (y_____ < height - 1) p2[width] = *p1;
			}
		}
		enumImageEnd;
	}


	//对区域进行标记
	template<typename T>
	int MarkRegion_MarkPixcel(T * bmp, int width, int height, int x, int y, int * region, int regionIndex)
	{
		int idx = width * y + x;
		if (!bmp[idx]) return 0;
		if (region[idx]) return 0;

		int area = 1;

		region[idx] = regionIndex;

		//switch (xmv::rand(4))
		//{
		//case 0:
		if (x > 0 && !region[idx - 1] && bmp[idx - 1])
			area += MarkRegion_MarkPixcel(bmp, width, height, x - 1, y, region, regionIndex) + 1;
		if (x < width - 1 && !region[idx + 1] && bmp[idx + 1])
			area += MarkRegion_MarkPixcel(bmp, width, height, x + 1, y, region, regionIndex) + 1;
		if (y > 0 && !region[idx - width] && bmp[idx - width])
			area += MarkRegion_MarkPixcel(bmp, width, height, x, y - 1, region, regionIndex) + 1;
		if (y < height - 1 && !region[idx + width] && bmp[idx + width])
			area += MarkRegion_MarkPixcel(bmp, width, height, x, y + 1, region, regionIndex) + 1;
		//break;
		//case 1:
		//	if (x < width - 1 && ! region[idx + 1] && bmp[idx + 1])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x + 1, y, region, regionIndex) + 1;
		//	if (x > 0 && ! region[idx - 1] && bmp[idx - 1])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x - 1, y, region, regionIndex) + 1;
		//	if (y < height - 1 && ! region[idx + width] && bmp[idx + width])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x, y + 1, region, regionIndex) + 1;
		//	if (y > 0 && ! region[idx - width] && bmp[idx - width])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x, y - 1, region, regionIndex) + 1;
		//	break;
		//case 2:
		//	if (y > 0 && ! region[idx - width] && bmp[idx - width])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x, y - 1, region, regionIndex) + 1;
		//	if (y < height - 1 && ! region[idx + width] && bmp[idx + width])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x, y + 1, region, regionIndex) + 1;
		//	if (x > 0 && ! region[idx - 1] && bmp[idx - 1])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x - 1, y, region, regionIndex) + 1;
		//	if (x < width - 1 && ! region[idx + 1] && bmp[idx + 1])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x + 1, y, region, regionIndex) + 1;
		//	break;
		//case 3:
		//	if (y < height - 1 && ! region[idx + width] && bmp[idx + width])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x, y + 1, region, regionIndex) + 1;
		//	if (y > 0 && ! region[idx - width] && bmp[idx - width])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x, y - 1, region, regionIndex) + 1;
		//	if (x < width - 1 && ! region[idx + 1] && bmp[idx + 1])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x + 1, y, region, regionIndex) + 1;
		//	if (x > 0 && ! region[idx - 1] && bmp[idx - 1])
		//		area += MarkRegion_MarkPixcel(bmp, width, height, x - 1, y, region, regionIndex) + 1;
		//	break;
		//}

		return area;
	}

	//对区域进行标记，区域序列号存储于region，区域大小存储于regionArea，返回Region的个数
	template<typename T>
	int MarkRegion(T * bmp, int * region, int * regionArea, int width, int height)
	{
		//删除一些区域
		T * ptr;
		int index = 1;
		int area = 0;
		int areaSum = 0;
		memset(region, 0, width * height * sizeof(int));
		enumImage(bmp, width, height, 0, width, 0, height, ptr)
		{
			if ((area = MarkRegion_MarkPixcel(bmp, width, height, x_____, y_____,
				region, index)) > 0)
			{
				regionArea[index++] = area;
				areaSum += area;
			}
		}
		enumImageEnd;
		regionArea[0] = width * height - areaSum;
		return index;
	}

	//删除指定的区域
	template<typename T>
	void DeleteRegion(T * bmp, int * region, int regionToDelete, int width, int height)
	{
		T * ptr;
		int  * r;
		enumArray2(bmp, region, 0, width * height, ptr, r)
		{
			if (*r == regionToDelete)
				*ptr = 0;
		}
		enumArrayEnd
	}


	//对图像重新采样，图像变为原来的1/n
	template<typename T>
	void ReSample_ReduceImageSize(T * bmp, int width, int height, T * newImage, int n, int newImageBufferSize, int & newImgWidth, int & newImgHeight)
	{
		newImgWidth = width / n;
		newImgHeight = height / n;

		T * pNewImage = newImage;
		T * p;
		int blockSize = n * n;
		int nx = 0;
		int ny = 0;
		for (int blocky = 0; ny < newImgHeight && blocky < height; ++ny, blocky += n)
		{
			nx = 0;
			for (int blockx = 0; blockx < width && nx < newImgWidth; ++nx, blockx += n, ++pNewImage)
			{
				T total = 0;
				enumImage(bmp, width, height, blockx, blockx + n, blocky, blocky + n, p)
				{
					total += *p;
				}
				enumImageEnd;
				*pNewImage = total /= blockSize;
			}
		}
	}
}