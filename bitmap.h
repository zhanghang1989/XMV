#pragma once
#include "rgb.h"
#include <fstream>
using namespace std;
#include "common.h"
#include "pixcel.h"

#include "string"
#include "matrix.h"
#include "array.h"
#include "vector.h"
#include "vectorH.h"


namespace xmv
{
	//template<typename T> struct BGR;
	//template<typename T> class Array;
	//template<typename T> class Matrix;
	//template<typename T> class Vector;
	//template<typename T> class VectorH;
	//create a Bitmap
	template<typename T> class Bitmap : public Matrix<BGR<T>>
	{
	public:
		Bitmap(string filename, T min = 0, T max = 255);			//Create a bitmap from file
		Bitmap(int width = 1, int hegiht = 1);
		Bitmap(Matrix<T> & gray);			//Create a bitmap from Matrix
#ifndef NORGB
		Bitmap(Matrix<RGB<T>> & m);			//Create a bitmap from Matrix
#endif
		Bitmap(Matrix<BGR<T>> & m);			//Create a bitmap from Matrix
		Bitmap(Matrix<T> & R, Matrix<T> & G, Matrix<T> & B);
		Bitmap(Vector<T> & gray, int width, int height);			//Create a bitmap from Vector
#ifndef NORGB
		Bitmap(Vector<RGB<T>> & m, int width, int height);			//Create a bitmap from Vector
#endif
		Bitmap(Vector<BGR<T>> & m, int width, int height);			//Create a bitmap from Vector
		Bitmap(Vector<T> & R, Vector<T> & G, Vector<T> & B, int width, int height);
		Bitmap(T * gray, int width, int height);
		Bitmap(T * R, T * G, T * B, int width, int height);
		Bitmap(BGR<T> * bgr, int width, int height);
#ifndef NORGB
		Bitmap(RGB<T> * rgb, int width, int height);
#endif

		T Min();
		T Max();

		Bitmap(Bitmap<T> & bmp);
		~Bitmap();

		int Width() { return Col(); }
		int Height() { return Row(); }

		void GetData(Matrix<BGR<T>> & m);
#ifndef NORGB
		void GetData(Matrix<RGB<T>> & m);
#endif
		void GetData(Matrix<T> & R, Matrix<T> & G, Matrix<T> & B);
		void GetData(Matrix<T> & gray);
		void GetData(Vector<BGR<T>> & m);
#ifndef NORGB
		void GetData(Vector<RGB<T>> & m);
#endif
		void GetData(Vector<T> & R, Vector<T> & G, Vector<T> & B);
		void GetData(Vector<T> & gray);
		void GetData(BGR<T> * buffer);
#ifndef NORGB
		void GetData(RGB<T> * buffer);
#endif
		void GetData(T * R, T * G, T * B);
		void GetData(T * gray);

		Bitmap<T> & operator = (Bitmap<T> & bmp);
		Bitmap<T> & operator = (Matrix<BGR<T>> & bmp);
		Bitmap<T> & operator = (Matrix<T> & bmp);
		Bitmap<T> & operator += (Bitmap<T> & bmp);
		BGR<T> * operator [] (int y);

		void Save(string filename, T min = 0, T max = 255);
		void FlipVertical();
		void FlipHorizontal();

		Bitmap<T> Cut(int xStart, int xSpan, int yStart, int ySpan);
		Matrix<T> Gray();

#ifdef CLR
	public:
		System::Windows::Forms::Form ^  ShowImage(T min = 0, T max = 255,
			System::Windows::Forms::IWin32Window ^ owner = nullptr,
			System::Windows::Forms::Form ^ form = nullptr
			);
		System::Drawing::Bitmap ^ ToBitmap(T min = 0, T max = 255);
		Bitmap(System::Drawing::Bitmap ^ bmp, T min = 0, T max = 255);
#endif
		//中值滤波
		Bitmap<T> MedianFilterSmoothing(int windowSize = 5, MedianFilterSmoothingWindowType windowType = MedianFilterSmoothingWindowType::Square);

		//获得HSV色彩空间
		void GetHSV(Matrix<T> & h, Matrix<T> & s, Matrix<T> & v);
	};

		template<typename T> Bitmap<T>::Bitmap(int width, int height)
			: Matrix(height, width)
		{
			memset(this->dataPtr, 0, this->memorySize);
		}

		template<typename T> Bitmap<T>::Bitmap(Matrix<T> & gray)
			: Matrix(gray.Row(), gray.Col())
		{
			GrayToRGB<T>(this->dataPtr, gray.Buffer(), gray.Col(), gray.Row());
		}

#ifndef NORGB
		template<typename T> Bitmap<T>::Bitmap(Matrix<RGB<T>> & m)		//Create a bitmap from Matrix
			: Matrix(m.Row(), m.Col())
		{
			RGBToRGB<T>(this->dataPtr, m.Buffer(), m.Col(), m.Row());
		}
#endif

		template<typename T> Bitmap<T>::Bitmap(Matrix<BGR<T>> & m)		//Create a bitmap from Matrix
			: Matrix(m.Row(), m.Col())
		{
			memcpy(this->dataPtr, m.Buffer(), m.MemorySize());
		}

		template<typename T> Bitmap<T>::Bitmap(Matrix<T> & R, Matrix<T> & G, Matrix<T> & B)
			: Matrix(R.Row(), R.Col())
		{
			RGBToRGB<T>(this->dataPtr, R.Buffer(), G.Buffer(), B.Buffer(), R.Col(), R.Row());
		}

		template<typename T> Bitmap<T>::Bitmap(Vector<T> & gray, int width, int height)
			: Matrix(height, width)
		{
			GrayToRGB<T>(this->dataPtr, gray.Buffer(), width, height);
		}

#ifndef NORGB
		template<typename T> Bitmap<T>::Bitmap(Vector<RGB<T>> & m, int width, int height)		//Create a bitmap from Vector
			: Matrix(height, width)
		{
			RGBToRGB<T>(this->dataPtr, m.Buffer(), width, height);
		}
#endif

		template<typename T> Bitmap<T>::Bitmap(Vector<BGR<T>> & m, int width, int height)		//Create a bitmap from Vector
			: Matrix(height, width)
		{
			memcopy(this->dataPtr, m.Buffer(), m.MemorySize());
		}

		template<typename T> Bitmap<T>::Bitmap(Vector<T> & R, Vector<T> & G, Vector<T> & B, int width, int height)
			: Matrix(height, width)
		{
			RGBToRGB<T>(this->dataPtr, R.Buffer(), G.Buffer(), B.Buffer(), width, height);
		}

		template<typename T> Bitmap<T>::Bitmap(T * gray, int width, int height)
			: Matrix(height, width)
		{
			GrayToRGB<T>(this->dataPtr, gray, width, height);
		}

		template<typename T> Bitmap<T>::Bitmap(T * R, T * G, T * B, int width, int height)
			: Matrix(height, width)
		{
			RGBToRGB<T>(this->dataPtr, R, G, B, width, height);
		}

		template<typename T> Bitmap<T>::Bitmap(BGR<T> * bgr, int width, int height)
			: Matrix(height, width)
		{
			memcpy(this->dataPtr, bgr, width * height * sizeof(BGR<T>));
		}

#ifndef NORGB
		template<typename T> Bitmap<T>::Bitmap(RGB<T> * rgb, int width, int height)
			: Matrix(height, width)
		{
			RGBToRGB<T>(this->dataPtr, rgb, width, height);
		}
#endif

		template<typename T> Bitmap<T>::Bitmap(Bitmap<T> & bmp)
			: Matrix(bmp.Height(), bmp.Width())
		{
			memcpy(this->dataPtr, bmp.Buffer(), bmp.MemorySize());
		}

		template<typename T> Bitmap<T>::~Bitmap()
		{
		}

		template<typename T> Bitmap<T> & Bitmap<T>::operator = (Bitmap<T> & bmp)
		{
			ReSize(bmp.Height(), bmp.Width());
			memcpy(this->dataPtr, bmp.Buffer(), bmp.MemorySize());
			return *this;
		}

		template<typename T> Bitmap<T> & Bitmap<T>::operator = (Matrix<BGR<T>> & bmp)
		{
			ReSize(bmp.Row(), bmp.Col());
			memcpy(this->dataPtr, bmp.Buffer(), bmp.MemorySize());
			return *this;
		}

		template<typename T> Bitmap<T> & Bitmap<T>::operator = (Matrix<T> & bmp)
		{
			ReSize(bmp.Row(), bmp.Col());

			T * p1; BGR<T> * p2;
			enumArray2(this->dataPtr, bmp.Buffer(), 0, this->size, p2, p1)
			{
				p2->B = p2->G = p2->R = *p1;
			}
			enumArrayEnd;

			return *this;
		}

		template<typename T> Bitmap<T> & Bitmap<T>::operator += (Bitmap<T> & bmp)
		{
			Bitmap<T> t = *this;
			BGR<T> * p1, *p2;
			enumArray2(this->Buffer(), bmp.Buffer(), 0, this->Size(), p1, p2)
			{
				*p1 += *p2;
			}
			enumArrayEnd;
			return *this;
		}

		template<typename T> BGR<T> * Bitmap<T>::operator[] (int y)
		{
			return this->dataPtr + y * col;
		}

		template<typename T> void Bitmap<T>::Save(string filename, T min = 0, T max = 255)
		{
			ofstream o(filename, std::ios_base::binary);

			string bmpHead = "BM";
			unsigned __int32 dataOffset = 54;
			unsigned __int32 headerSize = 40;
			__int32 width = Width();
			__int32 height = Height();
			unsigned __int16 bitsPerPixcel = 24;
			unsigned __int32 lineSize = (bitsPerPixcel * width / 32) * 4 + ((bitsPerPixcel * width % 32 != 0) ? 4 : 0);

			unsigned __int32 dataSize = lineSize * height;
			unsigned __int32 fileSize = headerSize + dataSize;

			WriteString(o, bmpHead);		//00 "BM"
			WriteUInt32(o, fileSize);		//02 4B 文件大小
			WriteZero(o, 4);				//06 4B 保留
			WriteUInt32(o, dataOffset);		//0A 4B 数据偏移地址

			WriteUInt32(o, headerSize);		//0E 4B 头结构大小
			WriteInt32(o, width);			//12 4B 宽度
			WriteInt32(o, height);			//16 4B 高度
			WriteUInt16(o, 1);				//1A 2B 色彩平面数，只能为1
			WriteUInt16(o, bitsPerPixcel);	//1C 2B 每个像素占用的比特数, 1, 4, 8, 16, 24, 32
			WriteUInt32(o, 0);				//1E 4B 压缩算法，不压缩为0
			WriteUInt32(o, 0);				//22 4B 图像大小
			WriteInt32(o, 1);				//26 4B 横向分辨率，像素/米
			WriteInt32(o, 1);				//2A 4B 纵向分辨率，像素/米
			WriteUInt32(o, 0);				//2E 4B 调色板颜色数，0表示为(2的色深次方)个
			WriteUInt32(o, 0);				//32 4B 重要颜色数量，0表示所有颜色都重要，一般不用本项

			//写入RGB数组
			//准备行数组
			byte * line = new byte[lineSize];
			memset(line, 0, lineSize);
			byte * ptr = line;
			BGR<T> * c = nullptr;
			BGR<T> * cy = this->dataPtr + width * (height - 1);
			T span = max - min;
			for (int y = height - 1; y >= 0; --y, cy -= width)
			{
				c = cy;
				ptr = line;
				for (int x = 0; x < width; ++x, ++c)
				{
					if (min == 0 && max == 255)
					{
						*ptr++ = c->G;
						*ptr++ = c->B;
						*ptr++ = c->R;
					}
					else
					{
						*ptr++ = (c->G - min) * 255 / span;
						*ptr++ = (c->B - min) * 255 / span;
						*ptr++ = (c->R - min) * 255 / span;
					}
				}
				//写入行
				WriteBytes(o, line, lineSize);
			}

			delete[]line;

			o.close();
		}

		template<typename T> Bitmap<T>::Bitmap(string filename, T min = 0, T max = 255)
		{
			ifstream i(filename, std::ios_base::binary);
			string bmpHead = "BM";

			string head = ReadString(i, 2);		//00 "BM"
			int fileSize = ReadUInt32(i);		//02 4B 文件大小
			ReadZero(i, 4);						//06 4B 保留
			int offSet = ReadUInt32(i);			//0A 4B 数据偏移地址
			int headSize = ReadUInt32(i);		//0E 4B 头结构大小
			int width = ReadInt32(i);			//12 4B 宽度
			int height = ReadInt32(i);			//16 4B 高度
			ReadUInt16(i);						//1A 2B 色彩平面数，只能为1
			int bitsPerPixcel = ReadUInt16(i);	//1C 2B 每个像素占用的比特数, 1, 4, 8, 16, 24, 32
			int cType = ReadUInt32(i);						//1E 4B 压缩算法，不压缩为0
			int imgSize = ReadUInt32(i);						//22 4B 图像大小
			int imgPX = ReadInt32(i);						//26 4B 横向分辨率，像素/米
			int imgPY = ReadInt32(i);						//2A 4B 纵向分辨率，像素/米
			int cpSize = ReadUInt32(i);			//2E 4B 调色板颜色数，0表示为(2的色深次方)个
			int icCount = ReadUInt32(i);		//32 4B 重要颜色数量，0表示所有颜色都重要，一般不用本项

			ReSize(height, width);
			unsigned __int32 lineSize = (bitsPerPixcel * width / 32) * 4 + ((bitsPerPixcel * width % 32 != 0) ? 4 : 0);

			if (offSet > 0x36)
			{
				byte * empty = new byte[offSet - 0x36];
				ReadBytes(i, empty, offSet - 0x36);
			}

			//写入RGB数组
			//准备行数组
			byte * line = new byte[lineSize];
			memset(line, 0, lineSize);

			BGR<T> * c = this->dataPtr;
			T span = max - min;

			if (bitsPerPixcel == 8)
			{

				for (int y = 0; y < height; ++y)
				{
					//读入
					ReadBytes(i, line, lineSize);
					byte * ptr = line;
					for (int x = 0; x < width; ++x, ++c)
					{
						if (min == 0 && max == 255)
						{
							c->G = c->B = c->R = *ptr++;
						}
						else
						{
							c->G = c->B = c->R = *ptr++ * span / 255 + min;
						}
					}
				}
			}


			if (bitsPerPixcel == 24)
			{

				for (int y = 0; y < height; ++y)
				{
					//读入
					ReadBytes(i, line, lineSize);
					byte * ptr = line;
					for (int x = 0; x < width; ++x, ++c)
					{
						if (min == 0 && max == 255)
						{
							c->G = *ptr++;
							c->B = *ptr++;
							c->R = *ptr++;
						}
						else
						{
							c->G = *ptr++ * span / 255 + min;
							c->B = *ptr++ * span / 255 + min;
							c->R = *ptr++ * span / 255 + min;
						}
					}
				}
			}

			delete[]line;

			i.close();

			this->FlipVertical();
		}

		template<typename T> void Bitmap<T>::GetData(Matrix<BGR<T>> & m)
		{
			m.ReSize(Height(), Width());
			RGBToRGB<T>(m.Buffer(), this->dataPtr, Width(), Height());
		}

#ifndef NORGB
		template<typename T> void Bitmap<T>::GetData(Matrix<RGB<T>> & m)
		{
			m.ReSize(Height(), Width());
			RGBToRGB<T>(m.Buffer(), this->dataPtr, Width(), Height());
		}
#endif

		template<typename T> void Bitmap<T>::GetData(Matrix<T> & R, Matrix<T> & G, Matrix<T> & B)
		{
			R.ReSize(Height(), Width());
			G.ReSize(Height(), Width());
			B.ReSize(Height(), Width());

			RGBToRGB(R.Buffer(), G.Buffer(), B.Buffer(), this->dataPtr, Width(), Height());
		}

		template<typename T> void Bitmap<T>::GetData(Matrix<T> & gray)
		{
			gray.ReSize(Height(), Width());
			RGBToGray<T>(gray.Buffer(), this->dataPtr, Width(), Height());
		}

		template<typename T> void Bitmap<T>::GetData(Vector<BGR<T>> & v)
		{
			v.ReSize(Width() * Height());
			RGBToRGB<T>(v.Buffer(), this->dataPtr, Width(), Height());
		}

#ifndef NORGB
		template<typename T> void Bitmap<T>::GetData(Vector<RGB<T>> & v)
		{
			v.ReSize(Width() * Height());
			RGBToRGB<T>(v.Buffer(), this->dataPtr, Width(), Height());
		}
#endif

		template<typename T> void Bitmap<T>::GetData(Vector<T> & R, Vector<T> & G, Vector<T> & B)
		{
			R.ReSize(Width() * Height());
			G.ReSize(Width() * Height());
			B.ReSize(Width() * Height());

			RGBToRGB(R.Buffer(), G.Buffer(), B.Buffer(), this->dataPtr, Width(), Height());
		}

		template<typename T> void Bitmap<T>::GetData(Vector<T> & gray)
		{
			gray.ReSize(Width() * Height());
			RGBToGray<T>(gray.Buffer(), this->dataPtr, Width(), Height());
		}

		template<typename T> void Bitmap<T>::GetData(BGR<T> * buffer)
		{
			memcpy(buffer, this->dataPtr, this->memorySize);
		}

#ifndef NORGB
		template<typename T> void Bitmap<T>::GetData(RGB<T> * buffer)
		{
			RGBToRGB<T>(buffer, this->dataPtr, Width(), Height());
		}
#endif

		template<typename T> void Bitmap<T>::GetData(T * R, T * G, T * B)
		{
			RGBToRGB<T>(R, G, B, this->dataPtr, Width(), Height());
		}

		template<typename T> void Bitmap<T>::GetData(T * gray)
		{
			RGBToGray<T>(gray, this->dataPtr, Width(), Height());
		}

#ifndef NORGB
		template<typename T>
		void RGBToRGB(T * destR, T * destG, T * destB, RGB<T> * sourceRgb, int width, int height)
		{
			RGB<T> * c;
			T * r, *g, *b;
			enumImage4(sourceRgb, destR, destG, destB, width, height, 0, width, 0, height, c, r, g, b)
			{
				*r = c->R;
				*g = c->G;
				*b = c->B;
			}
			enumImageEnd;
		}
#endif

		template<typename T>
		void RGBToRGB(T * destR, T * destG, T * destB, BGR<T> * sourceBgr, int width, int height)
		{
			BGR<T> * c;
			T * r, *g, *b;
			enumImage4(sourceBgr, destR, destG, destB, width, height, 0, width, 0, height, c, r, g, b)
			{
				*r = c->R;
				*g = c->G;
				*b = c->B;
			}
			enumImageEnd;
		}

#ifndef NORGB
		template<typename T>
		void RGBToRGB(RGB<T> * destRgb, T * sourceR, T * sourceG, T * sourceB, int width, int height)
		{
			RGB<T> * c;
			T * r, *g, *b;
			enumImage4(destRgb, sourceR, sourceG, sourceB, width, height, 0, width, 0, height, c, r, g, b)
			{
				c->R = *r;
				c->G = *g;
				c->B = *b;
			}
			enumImageEnd;
		}
#endif

		template<typename T>
		void RGBToRGB(BGR<T> * destBgr, T * sourceR, T * sourceG, T * sourceB, int width, int height)
		{
			BGR<T> * c;
			T * r, *g, *b;
			enumImage4(destBgr, sourceR, sourceG, sourceB, width, height, 0, width, 0, height, c, r, g, b)
			{
				c->R = *r;
				c->G = *g;
				c->B = *b;
			}
			enumImageEnd;
		}

#ifndef NORGB
		template<typename T>
		void RGBToRGB(BGR<T> * destBgr, RGB<T> * sourceRgb, int width, int height)
		{
			BGR<T> * bgr;
			RGB<T> * rgb;
			enumImage2(destBgr, sourceRgb, width, height, 0, width, 0, height,
				bgr, rgb)
			{
				bgr->R = rgb->R;
				bgr->G = rgb->G;
				bgr->B = rgb->B;
			}
			enumImageEnd;
		}

		template<typename T>
		void RGBToRGB(RGB<T> * destRgb, BGR<T> * sourceBgr, int width, int height)
		{
			BGR<T> * bgr;
			RGB<T> * rgb;
			enumImage2(destRgb, sourceBgr, width, height, 0, width, 0, height,
				rgb, bgr)
			{
				bgr->R = rgb->R;
				bgr->G = rgb->G;
				bgr->B = rgb->B;
			}
			enumImageEnd;
		}

		template<typename T>
		void RGBToGray(T * destGrayImg, RGB<T> * sourceRgbImg, int width, int height)
		{
			RGB<T> * rgb;
			T * t;
			enumImage2(sourceRgbImg, destGrayImg, width, height, 0, width, 0, height, rgb, t)
				*t = rgb->R * 0.299 + rgb->G * 0.587 + rgb->B * 0.114;
			enumImageEnd;
		}

		template<typename T>
		void RGBToGray(T * destGrayImg, BGR<T> * sourceBgrImg, int width, int height)
		{
			BGR<T> * rgb;
			T * t;
			enumImage2(sourceBgrImg, destGrayImg, width, height, 0, width, 0, height, rgb, t)
				*t = rgb->R * 0.299 + rgb->G * 0.587 + rgb->B * 0.114;
			enumImageEnd;
		}

		template<typename T>
		void GrayToRGB(RGB<T> * destRgbImg, T * sourceGrayImg, int width, int height)
		{
			RGB<T> * rgb;
			T * t;
			enumImage2(destRgbImg, sourceGrayImg, width, height, 0, width, 0, height, rgb, t)
				rgb->R = rgb->B = rgb->G = *t;
			enumImageEnd;
		}


#endif	
		template<typename T> void GrayToRGB(BGR<T> * destBgrImg, T * souceGrayImg, int width, int height)
		{
			BGR<T> * rgb;
			T * t;
			enumImage2(destBgrImg, souceGrayImg, width, height, 0, width, 0, height, rgb, t)
				rgb->R = rgb->B = rgb->G = *t;
			enumImageEnd;
		}


		template<typename T> Bitmap<T> Bitmap<T>::Cut(int xStart, int xSpan, int yStart, int ySpan)
		{
			Bitmap<T> m(xSpan, ySpan);
			BGR<T> * ptrD = m.Buffer();
			BGR<T> * ptrS = this->dataPtr + col * yStart + xStart;
			int copyLen = xSpan * sizeof(BGR<T>);
			for (int i = 0; i < ySpan; ++i, ptrD += xSpan, ptrS += this->col)
				memcpy(ptrD, ptrS, copyLen);

			return m;
		}	 //获得矩阵的一部分

		template<typename T> T Bitmap<T>::Min()
		{
			T m = this->dataPtr->Min();
			for (BGR<T> * ptr = dataPtr; ptr < dataEndPtr; ++ptr)
			{
				T v = ptr->Min();
				if (v < m) m = v;
			}

			return m;
		}

		template<typename T> T Bitmap<T>::Max()
		{
			T m = this->dataPtr->Max();
			for (BGR<T> * ptr = dataPtr; ptr < dataEndPtr; ++ptr)
			{
				T v = ptr->Max();
				if (v > m) m = v;
			}

			return m;
		}


#ifdef CLR
		template <typename T>
		System::Windows::Forms::Form ^ Bitmap<T>::ShowImage(T min, T max,
			System::Windows::Forms::IWin32Window ^ owner,
			System::Windows::Forms::Form ^ form)
		{
			//System::String ^ tempPath = System::Windows::Forms::Application::UserAppDataPath;
			//tempPath += "\\" + System::Guid::NewGuid().ToString() + ".bmp";

			//this->Save(ToString(tempPath), min, max);

			System::Drawing::Bitmap ^ bmp = ToBitmap(min, max);
			if (form == nullptr || form->IsDisposed)
				form = gcnew System::Windows::Forms::Form();
			form->ClientSize = bmp->Size;

			System::Windows::Forms::PictureBox ^ box;
			if (form->Controls->Count == 0)
			{
				box = gcnew System::Windows::Forms::PictureBox();
				form->Controls->Add(box);
			}
			else
				box = (System::Windows::Forms::PictureBox ^)form->Controls[0];

			box->SizeMode = System::Windows::Forms::PictureBoxSizeMode::AutoSize;
			box->Dock = System::Windows::Forms::DockStyle::Fill;
			box->Image = gcnew System::Drawing::Bitmap(bmp);

			form->AutoScroll = true;

			form->Show();
			form->Refresh();
			System::Windows::Forms::Application::DoEvents();

			//System::IO::File::Delete(tempPath);

			//f->Show(nullptr);
			return form;
		}

		template <typename T>
		System::Drawing::Bitmap ^ Bitmap<T>::ToBitmap(T min, T max)
		{
			int width = Width();
			int height = Height();
			System::Drawing::Bitmap ^ bmp = gcnew System::Drawing::Bitmap(width, height);
			System::Drawing::Imaging::BitmapData ^ data = bmp->LockBits(System::Drawing::Rectangle(0, 0, width, height),
				System::Drawing::Imaging::ImageLockMode::WriteOnly, System::Drawing::Imaging::PixelFormat::Format24bppRgb);

			byte * buff;
			if (sizeof(buff) == 8)
				buff = (byte *)(data->Scan0.ToInt64());
			else
				buff = (byte *)(data->Scan0.ToInt32());

			int cw = data->Stride;

			BGR<T> * c = this->dataPtr;
			T span = max - min;
			for (int y = 0; y < height; ++y)
			{
				byte * ptr = buff + y * cw;
				for (int x = 0; x < width; ++x, ++c)
				{
					if (min == 0 && max == 255)
					{
						*ptr++ = c->G;
						*ptr++ = c->B;
						*ptr++ = c->R;
					}
					else
					{
						*ptr++ = (c->G - min) * 255 / span;
						*ptr++ = (c->B - min) * 255 / span;
						*ptr++ = (c->R - min) * 255 / span;
					}
				}
				//写入行
			}

			bmp->UnlockBits(data);

			return bmp;
		}

		template<typename T> Bitmap<T>::Bitmap(System::Drawing::Bitmap ^ bmp, T min = 0, T max = 255)
		{
			int width = bmp->Width;
			int height = bmp->Height;

			Matrix::ReSize(height, width);

			System::Drawing::Imaging::BitmapData ^ data = bmp->LockBits(System::Drawing::Rectangle(0, 0, width, height),
				System::Drawing::Imaging::ImageLockMode::ReadOnly, System::Drawing::Imaging::PixelFormat::Format24bppRgb);
#ifdef X64
			byte * buff = (byte *)(data->Scan0.ToInt64());
#else
			byte * buff = (byte *)(data->Scan0.ToInt32());
#endif
			int cw = data->Stride;

			BGR<T> * c = this->dataPtr;
			T span = max - min;
			for (int y = 0; y < height; ++y)
			{
				byte * ptr = buff + y * cw;
				for (int x = 0; x < width; ++x, ++c)
				{
					if (min == 0 && max == 255)
					{
						c->G = *ptr++;
						c->B = *ptr++;
						c->R = *ptr++;
					}
					else
					{
						c->G = *ptr++ * span / 255 + min;
						c->B = *ptr++ * span / 255 + min;
						c->R = *ptr++ * span / 255 + min;
					}
				}
				//写入行
			}

			bmp->UnlockBits(data);
		}
#endif

		//图像垂直翻转
		template <typename T> void Bitmap<T>::FlipVertical()
		{
			int h = Height();
			for (int s = 0, e = h - 1; s < e; ++s, --e)
				TransformRowSwap(s, e);
		}

		//图像水平翻转
		template <typename T> void Bitmap<T>::FlipHorizontal()
		{
			int h = Width();
			for (int s = 0, e = h - 1; s < e; ++s, --e)
				TransformColSwap(s, e);
		}

		template <typename T> Matrix<T> Bitmap<T>::Gray()
		{
			Matrix<T> m;
			this->GetData(m);
			return m;
		}

		template <typename T> Bitmap<T> Bitmap<T>::MedianFilterSmoothing(int windowSize, MedianFilterSmoothingWindowType windowType)
		{
			Matrix<T> R, G, B;
			this->GetData(R, G, B);

			R = R.MedianFilterSmoothing(windowSize, windowType);
			G = G.MedianFilterSmoothing(windowSize, windowType);
			B = B.MedianFilterSmoothing(windowSize, windowType);

			return Bitmap<T>(R, G, B);
		}

		template<typename T> void Bitmap<T>::GetHSV(Matrix<T> & h, Matrix<T> & s, Matrix<T> & v)
		{
			int row = Row();
			int col = Col();
			h.ReSize(row, col);
			s.ReSize(row, col);
			v.ReSize(row, col);

			BGR<T> * p;
			T * ph, *ps, *pv;
			enumArray4(this->Buffer(), h.Buffer(), s.Buffer(), v.Buffer(), 0, this->Size(),
				p, ph, ps, pv)
			{
				T max = xmv::Max(xmv::Max(p->R, p->G), p->B);
				T min = xmv::Min(xmv::Min(p->R, p->G), p->B);
				if (max == min)
					*ph = 0;
				else
				{
					if (p->R == max) *ph = (p->G - p->B) / (max - min);
					if (p->G == max) *ph = 2 + (p->B - p->R) / (max - min);
					if (p->B == max) *ph = 4 + (p->R - p->G) / (max - min);

					*ph *= 60;

					if (*ph < 0) *ph = *ph + 360;
					*ph /= 360;
				}

				*pv = max;
				*ps = (max - min) / max;
			}
			enumArrayEnd;
		}
}