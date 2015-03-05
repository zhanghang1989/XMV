#pragma once
#include <convert.h>
#include <bitmap.h>
#include <matrix.h>
#include <kernel.h>

namespace xmv
{
	class HOG
	{
	public:
		HOG(int bsize = 16, int bnum = 9);
		HOG(int bsize, int bnum, double x, double y, double r);
		HOG(int bsize, int bnum, double x, double y, double w, double h);

		HOG(string filename, int bsize = 16, int bnum = 9);
		HOG(string filename, int bsize, int bnum, double x, double y, double r);
		HOG(string filename, bool default_cicle, int bsize = 16, int bnum = 9);
		HOG(string filename, int bsize, int bnum, double x, double y, double w, double h);

		HOG(Bitmap<byte> img, int bsize = 16, int bnum = 9);
		HOG(Bitmap<byte> img, bool default_cicle, int bsize = 16, int bnum = 9);
		HOG(Bitmap<byte> img, int bsize, int bnum, double x, double y, double r);
		HOG(Bitmap<byte> img, int bsize, int bnum, double x, double y, double w, double h);
		~HOG();

	private:
		// intial parameters
		int blksize;
		int binnum;
		int width;                  
		int height;

		int m_feavectorsize;
		int blknum;
		int blkvectorsize;

		// Rgb and Gray image
		Bitmap<byte> rgb_img;
		Bitmap<byte> gray_img;

		// window of weight
		Matrix<double> m_gaussian;

		// Gradient
		Matrix<BGR<double>> m_gradient[2];
		Matrix<double> n_gradient;
		Matrix<double> m_angle;

		// Hog features
		Matrix<double> hogfeature;
		vector<Point<int>> featurecenter;

		// Gamma table 
		byte GAMMATABLE[256];
		const float GAMMA = 0.5;
		const float NORMALOFFSET = 0.5;
		const int OVERLAP = 2;

		// ROI
		bool UseROI;
		// circle = true, rectangle =false
		bool CircleROI;
		double roixs, roiys, roiwidth, roiheight, roiradius;

	public:

		// // intialization
		void Init(string);
		void Init(string filename, double x, double y, double r);
		void Init(string filename, double x, double y, double w, double h);
		void Init(Bitmap<byte> img);
		void Init(Bitmap<byte> img, double x, double y, double r);
		void Init(Bitmap<byte> img, double x, double y, double w, double h);
		

		// Solving features
		void BuildGammaTable(float);
		void GammaNormAndCorrection();
		void SolveGradient();
		Matrix<double> HogFeature();
		void SaveFeature2Img(std::string);
		void SetROI(double, double, double, double);
		void SetROI(double, double, double);

		// Test and Debug
		template<typename T>
		void MarkCross(BGR<T>* buffer, int width, int height, int x, int y, int size = 3, BGR<T> color = BGR<T>(0,0,255));

		Matrix<double> GetFeature()
		{
			return hogfeature;
		}
		int Get_FeaVectorSize()
		{
			return m_feavectorsize;
		}
		int Getblknum()
		{
			return blknum;
		}
		int Getblkvectorsize()
		{
			return blkvectorsize;
		}
		vector<Point<int>> GetFeaturePoints()
		{
			return featurecenter;
		}

		void test()
		{
			gray_img.ShowImage();
			Bitmap<double> grad_x(m_gradient[0]);
			Bitmap<double> grad_y(m_gradient[1]);
			Bitmap<double> scale(n_gradient);
			Bitmap<double> angle(m_angle);
			grad_x.ShowImage();
			grad_y.ShowImage();
			cout << scale.Max();
			scale.ShowImage();
			angle.ShowImage();
			Bitmap<double> outimg(hogfeature);
			outimg = outimg * 255;
			outimg.ShowImage();
			Array<int> hist = m_angle.histgram(9);
			hist.WriteToTextFile("hist.txt");
			ofstream out("hogfeatures.txt");
			out << hogfeature;
			out.close();
			system("pause");
		}
	};

	HOG::HOG(int bsize, int bnum)
		:blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
	}

	HOG::HOG(string filename, int bsize, int bnum)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		UseROI = false;
		Init(filename);
	}

	HOG::HOG(Bitmap<byte> img, int bsize, int bnum)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		UseROI = false;
		Init(img);
	}

	HOG::HOG(int bsize, int bnum, double x, double y, double r)
		:blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		SetROI(x, y, r);
	}

	HOG::HOG(int bsize, int bnum, double x, double y, double w, double h)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		SetROI(x, y, w, h);
	}

	HOG::HOG(string filename, int bsize, int bnum, double x, double y, double r)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		Init(filename,x, y, r);
	}

	HOG::HOG(string filename, int bsize, int bnum, double x, double y, double w, double h)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		Init(filename, x, y, w, h);
	}

	HOG::HOG(Bitmap<byte> img, int bsize, int bnum, double x, double y, double r)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		Init(img, x, y, r);
	}

	HOG::HOG(Bitmap<byte> img, int bsize, int bnum, double x, double y, double w, double h)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		Init(img, x, y, w, h);
	}

	HOG::HOG(string filename, bool default_cicle, int bsize, int bnum)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		Bitmap<byte> img(filename);
		int x = img.Width() / 2;
		int y = img.Height() / 2;
		int r = y - 2;
		Init(filename, x, y, r);
	}

	HOG::HOG(Bitmap<byte> img, bool default_cicle, int bsize, int bnum)
		: blksize(bsize), binnum(bnum), blknum(0), blkvectorsize(0)
	{
		int x = img.Width() / 2;
		int y = img.Height() / 2;
		int r = y - 1;
		Init(img, x, y, r);
	}

	HOG::~HOG()
	{
	}

	void HOG::Init(string filename)
	{
		Bitmap<byte> img(filename);
		rgb_img = img;
		gray_img = img.Gray();

		BuildGammaTable(GAMMA);
		width = gray_img.Width();
		height = gray_img.Height();
		m_gaussian = CreateGaussianFilterKernel<double>(blksize);
		HogFeature();
	}

	void HOG::Init(Bitmap<byte> img)
	{
		rgb_img = img;
		gray_img = img.Gray();

		BuildGammaTable(GAMMA);
		width = gray_img.Width();
		height = gray_img.Height();
		m_gaussian = CreateGaussianFilterKernel<double>(blksize);
		HogFeature();
	}

	void HOG::Init(string filename, double x, double y, double r)
	{
		SetROI(x, y, r);
		Init(filename);
	}

	void HOG::Init(string filename, double x, double y, double w, double h)
	{
		SetROI(x, y, w, h);
		Init(filename);
	}

	void HOG::Init(Bitmap<byte> img, double x, double y, double r)
	{
		SetROI(x, y, r);
		Init(img);
	}

	void HOG::Init(Bitmap<byte> img, double x, double y, double w, double h)
	{
		SetROI(x, y, w, h);
		Init(img);
	}

	void HOG::BuildGammaTable(float pGamma)
	{
		float f;
		for (int i = 0; i < 256; i++)
		{
			f = (i + NORMALOFFSET) / 256;
			f = (float) pow(f, pGamma);
			GAMMATABLE[i] = (byte) (f * 256 - NORMALOFFSET);
		}
	}

	void HOG::GammaNormAndCorrection()
	{
		int img_size = gray_img.Width() * gray_img.Height();
		BGR<byte> *img_ptr = gray_img.Buffer();
		for (int i = 0; i < img_size; i++)
		{
			*(img_ptr + i) = GAMMATABLE[(img_ptr + i)->B];
		}
	}

	void HOG::SolveGradient()
	{
		m_gradient[0].ReSize(gray_img.Height(), gray_img.Width());
		m_gradient[1].ReSize(gray_img.Height(), gray_img.Width());
		// Matrix<double> A

		n_gradient.ReSize(gray_img.Height(), gray_img.Width());
		m_angle.ReSize(gray_img.Height(), gray_img.Width());

		// x-direction gradient
		SampleLoader<BGR<double>> sload1;
		Vector<BGR<double>> data(gray_img.Height(), 0);
		for (int i = 1; i < gray_img.Width(); i++)
		{
			sload1.AddSample(gray_img.GetCol(i));
		}
		sload1.AddSample(data);
		Matrix<BGR<double>> left = sload1.GetSampleMatrix();

		SampleLoader<BGR<double>> sload2;
		sload2.AddSample(data);
		for (int i = 0; i < gray_img.Width() - 1; i++)
		{
			sload2.AddSample(gray_img.GetCol(i));
		}
		Matrix<BGR<double>> right = sload2.GetSampleMatrix();

		m_gradient[0] = left - right;

		// y-direction gradient
		SampleLoader<BGR<double>> sload3;
		Vector<BGR<double>> data1(gray_img.Width(), 0);
		for (int i = 1; i < gray_img.Height(); i++)
		{
			sload3.AddSample(gray_img.GetRow(i).Transpose());
		}
		sload3.AddSample(data1);
		Matrix<BGR<double>> top = sload3.GetSampleMatrix().Transpose();

		SampleLoader<BGR<double>> sload4;
		sload4.AddSample(data1);
		for (int i = 0; i < gray_img.Height() - 1; i++)
		{
			sload4.AddSample(gray_img.GetRow(i).Transpose());
		}
		Matrix<BGR<double>> down = sload4.GetSampleMatrix().Transpose();

		m_gradient[1] = top - down;
		// cout << "m_grad = " << m_gradient[0].Row() << '\t' << m_gradient[0].Col() << endl;
		// cout << "m_grad = " << m_gradient[1].Row() << '\t' << m_gradient[1].Col() << endl;
		
		// calculate the magnitude and angle

		BGR<double> *ptr1 = m_gradient[0].Buffer();
		BGR<double> *ptr2 = m_gradient[1].Buffer();
		double * ptr3 = n_gradient.Buffer();
		double *ptr_angle = m_angle.Buffer();

		for (int i = 0; i < gray_img.Size(); i++)
		{
			*(ptr3 + i) = sqrt(1.0 * (ptr1 + i)->B * (ptr1 + i)->B + (ptr2 + i)->B * (ptr2 + i)->B);
			if ((ptr1 + i)->B == 0)
				*(ptr_angle + i) = EPS;
			else
				*(ptr_angle + i) = (atan((ptr2 + i)->B/ (ptr1 + i)->B) + PI/2)/ PI * 180;
		}
	}

	Matrix<double> HOG::HogFeature()
	{
		GammaNormAndCorrection();
		SolveGradient();

		// store each block histogram
		double cellsize = blksize / 2;

		int step = blksize / OVERLAP;
		int xEnd = width - blksize + 1;
		int yEnd = height - blksize + 1;

		int xblknum = (width - blksize) / step + 1;
		int yblknum = (height - blksize) / step + 1;

		int blkbuffersize = 4 * binnum;
		m_feavectorsize = xblknum * yblknum * blkbuffersize;
		blknum = xblknum * yblknum;

		bool interpolation = true;
		int m, n;
		int indexX, indexY, indexZ, indexXN, indexYN, indexZN;
		double weight;
		double gradient_orient, gradient_scale;
		int halfangle = 90 / binnum;
		int fullangle = 180 / binnum;
		double  x1, y1, z1, x2, y2;
		int i, j, bi, bj, x, y;
		double ix, jy, gz, ix1, jy1;

		hogfeature.ReSize(blknum, blkbuffersize);
		VectorH<double> blockHist(blkbuffersize);
		Array<double> sume(binnum);

		int blocknum = 0;

		for (i = 0; i < yEnd; i += step)
		{
			for (j = 0; j < xEnd; j += step)
			{
				if (UseROI && CircleROI)
				{
					int cord [] = { i + cellsize - roiys, j + cellsize - roixs };
					Array<double> cord_ar(2, cord, 2);
					if (cord_ar.Norm2() >= (roiradius - 0.7070 * cellsize))
						continue;
				}
				else if (UseROI && !CircleROI)
				{
					if (i < roiys || j < roixs || i >(roiys + roiheight) || j >(roixs + roiwidth))
						continue;
				}
				featurecenter.push_back(Point<int>(j, i));
				BGR<byte> color(0, 0, 255);
				MarkCross(gray_img.Buffer(), gray_img.Width(), gray_img.Height(), j + cellsize, i + cellsize);// , 5, color);

				blockHist = 0;
				for (bi = 0; bi < blksize; bi++)
				{
					for (bj = 0; bj < blksize; bj++)
					{
						interpolation = true;
						m = i + bi;
						n = j + bj;

						weight = m_gaussian.GetValue(bi, bj);
						gradient_scale = n_gradient.GetValue(m, n);
						gradient_orient = m_angle.GetValue(m, n);

						// find the cell index which nearest to current point
						// indexX for 1-d, indexY for 2-d, indexZ for 3-d(theta)
						indexX = floor((bi + 1 + cellsize) / blksize);
						indexY = floor((bj + 1 + cellsize) / blksize);
						indexZ = (int) (gradient_orient / fullangle) % binnum;

						x1 = indexX*cellsize + cellsize / 2 - 1;
						y1 = indexY*cellsize + cellsize / 2 - 1;
						z1 = indexZ*fullangle + halfangle;


						// find symmetric cell index
						indexXN = 1 - indexX;
						indexYN = 1 - indexY;

						x2 = indexXN*cellsize + cellsize / 2 - 1;
						y2 = indexYN*cellsize + cellsize / 2 - 1;

						if (indexZ != 0 && indexZ != binnum - 1)
						{
							indexZN = (gradient_orient - indexZ*fullangle) >= halfangle ? indexZ + 1 : indexZ - 1;
						}
						else
						{
							interpolation = false;
						}
						// interpolation = false;

						if (interpolation)
						{
							ix = abs(bi - x1);
							jy = abs(bj - y1);
							gz = abs(gradient_orient - z1);

							ix1 = abs(bi - x2);
							jy1 = abs(bj - y2);
							if (ix > cellsize || jy > cellsize || ix1 > cellsize || jy1 > cellsize)
								interpolation = false;
						}
						// change from dalaa's paper for improving the code efficiency
						// find the cell which has 8-neighbour, then we can do trilinear interpolation
						if (interpolation)
						{

							// trilinear interpolation
							int cnt1 = indexX * 2 * binnum + indexY * binnum + indexZ;
							int cnt2 = indexX * 2 * binnum + indexY * binnum + indexZN;

							blockHist[cnt1] += gradient_scale*weight
								*(1.00 - ix / cellsize)*(1.00 - jy / cellsize)*(1.00 - gz / fullangle);

							blockHist[cnt2] += gradient_scale*weight
								* (1.00 - ix / cellsize)*(1.00 - jy / cellsize)*(gz / fullangle);

							blockHist[cnt1] += gradient_scale*weight
								* (1.00 - ix / cellsize)*(jy / cellsize)*(1.00 - gz / fullangle);


							blockHist[cnt2] += gradient_scale*weight
								* (ix / cellsize)*(1.00 - jy / cellsize)*(1.00 - gz / fullangle);

							blockHist[cnt1] += gradient_scale*weight
								* (ix / cellsize)*(jy / cellsize)*(1.00 - gz / fullangle);

							blockHist[cnt2] += gradient_scale*weight
								* (1.00 - ix / cellsize)*(jy / cellsize)*(gz / fullangle);

							blockHist[cnt1] += gradient_scale*weight
								* (ix / cellsize)*(1.00 - jy / cellsize)*(gz / fullangle);

							blockHist[cnt2] += gradient_scale*weight
								* (ix / cellsize)*(jy / cellsize)*(gz / fullangle);
						}
						else
						{
							int cnt = indexX * 2 * binnum + indexY * binnum + indexZ;
							blockHist[cnt] += gradient_scale*weight;
						}
					}
				}
				// normalize feature 
				sume = 0;

				for (x = 0; x < binnum; x++)
				{
					for (y = 0; y < 4; y++)
					{
						sume[x] += blockHist[y*binnum + x] * blockHist[y*binnum + x];
					}
					sume[x] = sqrt(sume[x] + EPS);
					for (y = 0; y < 4; y++)
					{
						blockHist[y*binnum + x] = blockHist[y*binnum + x] / sume[x];
					}
				}

				hogfeature.SetRow(blocknum, blockHist);
				blocknum++;
			}
		}
		blknum = blocknum;
		blkvectorsize = blkbuffersize;

		if (UseROI)
		{
			hogfeature = hogfeature.Cut(0, blknum, 0, blkvectorsize);
		}


		return hogfeature;
	}

	void HOG::SaveFeature2Img(std::string filename)
	{
		Bitmap<double> outimg(hogfeature);
		outimg = outimg * 255;
		outimg.Save(filename);
	}

	void HOG::SetROI(double x, double y, double radius)
	{
		UseROI = true;
		CircleROI = true;
		roixs = x;
		roiys = y;
		roiradius = radius;
	}

	void HOG::SetROI(double x, double y, double w, double h)
	{
		UseROI = true;
		CircleROI = false;
		roixs = x;
		roiys = y;
		roiwidth = w;
		roiheight = h;
	}

	template <typename T>
	void HOG::MarkCross(BGR<T>* ptr, int width, int height, int x, int y, int size, BGR<T> color)
	{
		for (int i = x - size; i < x + size + 1; i++)
		{
			*(ptr + y * width + i) = color;
		}
		for (int i = y - size; i < y + size + 1; i++)
		{
			*(ptr + i * width + x) = color;
		}
	}
}
