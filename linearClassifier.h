#pragma once

#include <fstream>
#include "vector.h"
#include "matrix.h"
#include <list>
#include <cstring>

using namespace std;

namespace xmv
{
	template<typename T> struct AdaBoostStepValue
	{
		T step;
		long long usedCount;
	};

	template<typename T> int AdaBoostStepValueCompare(AdaBoostStepValue<T> & x, AdaBoostStepValue<T> & y)
	{
		T v1 = x.step * x.usedCount;
		T v2 = y.step * y.usedCount;
		if (v1 == v2) return 0;
		else if (v1 > v2) return -1;
		else return 1;
	}

	template<typename T> class LinearClassifier
	{
	private:
		list<T*> samples;
		list<T> classValues;
		Vector<T> parameters;
		Matrix<T> abParameters;
		Vector<T> abPowers;

		double move;

		bool mul;
		int parameterCount;
		int power;
	public:
		LinearClassifier(int parameterCount, int power = 1, bool mul = false)
			: parameters(parameterCount * power + (mul ? parameterCount : 0) + 1, (T)0)
		{
			this->power = power;
			this->parameterCount = parameterCount;
			this->mul = mul;
			this->move = 0;
		}

		LinearClassifier(int power = 1, bool mul = false)
		{
			this->power = power;
			this->parameterCount = -1;
			this->mul = mul;
			this->move = 0;
		}

		void Move(double m)
		{
			move = m;
		}

		int ParameterSize()
		{
			return parameterCount * power + (mul ? parameterCount : 0) + 1;
		}

		~LinearClassifier()
		{
			for (list<T*>::iterator p = samples.begin(); p != samples.end(); ++p)
				DELETEARR(*p);
		}

	public:
		void AddSample(Vector<T> & sample, T classValue)
		{
			if (this->parameterCount <= 0)
			{
				parameterCount = sample.Size();
				parameters = Vector<T>(parameterCount * power + (mul ? parameterCount : 0) + 1, (T)0);
			}

			T* arr = NEW("XMV.LinearClassifier.AddSample: arr") T[this->parameterCount];
			::memcpy(arr, sample.Buffer(), parameterCount * sizeof(T));
			samples.push_back(arr);
			classValues.push_back(classValue);
		}

	protected:
		void GetSampleAndResultMatrix(Matrix<T> & m, Vector<T> & v)
		{
			m.ReSize(samples.size(), ParameterSize());
			v.ReSize(samples.size());

			int rowCount = samples.size();

			T * mPtr = m.Buffer();
			T * vPtr = v.Buffer();

			for (list<T *>::iterator samplePtr = samples.begin(); samplePtr != samples.end(); ++samplePtr)
			{
				T * valuePtr;
				T first = **samplePtr;
				T * lineEnd = *samplePtr + parameterCount;
				enumArray(*samplePtr, 0, parameterCount, valuePtr)
				{
					T v = 1;
					for (int k = 1; k <= power; ++k)
						*mPtr++ = v *= *valuePtr;

					if (mul)
					{
						if (valuePtr + 1 < lineEnd)
							*mPtr++ = *valuePtr * valuePtr[1];
						else
							*mPtr++ = *valuePtr * first;
					}
				} enumArrayEnd;

				*mPtr++ = 1;
			};

			for (list<T>::iterator p = classValues.begin(); p != classValues.end(); ++p)
			{
				*vPtr++ = *p;
			};

		}

	public:
		//最小二乘法学习
		void LeastSquaresMethodLearn()
		{
			Matrix<T> m;
			Vector<T> v;

			GetSampleAndResultMatrix(m, v);

			this->parameters = SolveApproximateLinearEquations(m, v);
		}

	private:
		//在不进行矩阵乘法的情况下进行快速判别梯度下降的优劣
		int FastTest(Matrix<T> & samples, Vector<T> & x, Vector<T> & sampleX, Vector<T> & b, int changingId, T changingValue, T & ERVal)
		{
			Vector<T> k = sampleX + samples.GetCol(changingId) * changingValue;
			T * pk, *pb;
			T err = 0;
			ERVal = 0;

			enumArray2(k.Buffer(), b.Buffer(), 0, k.Size(), pk, pb)
			{
				if (*pk < 0 && *pb >= 0 || *pk >= 0 && *pb < 0)
				{
					++err;
					ERVal += Abs<T>(*pb);
				}
			}
			enumArrayEnd;

			return err;
		}

	private:
		//测试学习误差
		int Test(Matrix<T> & samples, Vector<T> & x, Vector<T> & sampleX, Vector<T> & b, T & ERVal, T errIncreacement = 0)
		{
			Vector<T> & k = (sampleX = samples * x);
			T * pk, *pb;
			T err = 0;
			ERVal = 0;

			enumArray2(k.Buffer(), b.Buffer(), 0, k.Size(), pk, pb)
			{
				if (*pk < 0 && *pb >= 0 || *pk >= 0 && *pb < 0)
				{
					++err;
					ERVal += Abs<T>(*pb);

					*pb *= (1 + errIncreacement);
					//if (*pb < 0) 
					//	*pb -= errIncreacement;
					//else 
					//	*pb += errIncreacement;
				}
			}
			enumArrayEnd;

			return err;
		}

	public:
		//梯度下降学习
		//stepStart: 起始步长，stepIncreasement：递增步长，stepMax：最大步长
		Vector<T> GradientdescentLearn(T stepStart = 0.001, T stepIncreasement = 0.001, T stepMax = 0.1)
		{
			Matrix<T> m;
			Vector<T> v;
			GetSampleAndResultMatrix(m, v);

			Vector<T> x = SolveApproximateLinearEquations(m, v);

			return GradientdescentLearn(m, x, v, stepStart, stepIncreasement, stepMax);
		}

	private:
		Matrix<T> elmA;
		Vector<T> elmB;
		Vector<T> elmX;


	public:
		//极限学习
		void ELMLearn(string path, int hiddenCount, T stepStart = 0.001, T stepIncreasement = 0.001, T stepMax = 0.1)
		{
			//创建矩阵H
			Matrix<T> m;
			Vector<T> v;
			GetSampleAndResultMatrix(m, v);
			int sampleCount = m.Row();
			int sampleSize = m.Col();

			Matrix<T> H(sampleCount, hiddenCount);

			elmA = RandomMatrix<T>(sampleSize, hiddenCount, -1, 1);
			//elmB = RandomVector<T>(hiddenCount, -1, 1);
			elmX = Vector<T>(hiddenCount, 0);
			SaveElmParameters(path);

			H = m * elmA;
			T * h = H.Buffer();
			for (int y = 0; y < sampleCount; ++y)
			{
				//T * b = elmB.Buffer();
				for (int x = 0; x < hiddenCount; ++x, /*++b,*/ ++h)
				{
					*h = tanh(*h);
				}
			}

			H = MergeMatrixHorizontal(H, Matrix<T>(sampleCount, 1, 1));
			SaveElmParameters(path);
			cout << "H Matrix OK!                  \r";

			Vector<T> x = SolveApproximateLinearEquations(H, v);

			elmX = GradientdescentLearn(H, x, v, stepStart, stepIncreasement, stepMax);
			SaveElmParameters(path);
		}

	public:
		void SaveElmParameters(string filename)
		{

			ofstream o(filename);
			o << elmA.Row() << endl;
			o << elmA.Col() << endl;
			o << format(8, 10) << elmA << endl;

			//o << elmB.Size() << endl;
			//o << format(8, 10) << elmB << endl;

			o << elmX.Size() << endl;
			o << format(8, 10) << elmX << endl;
			o.close();
		}

	private:
		Vector<T> GradientdescentLearn(Matrix<T> m, Vector<T> initX, Vector<T>v, AdaBoostStepValue<T> * steps, int stepCount, int maxN = -1, int nDirection = 1)
		{
			int stepIndex = 0;
			T step = steps[stepIndex].step;

			Vector<T> x = initX;
			//SolveApproximateLinearEquations(m, v);
			bool lastIsOK = false;

			int n = maxN < 0 ? x.Size() : maxN;
			T ErV = -1;
			for (;;)
			{
				Vector<T> mx;
				T OutErV;
				int Er = Test(m, x, mx, v, OutErV);
				if (ErV == -1 || OutErV <= ErV)
					ErV = OutErV;
				cout << Er << " / " << m.Row() << " | ERV:" << ErV << ", Step:" << step << "/" << stepIndex << "        \r";

				int cc = 0;

				if (nDirection == -1)
				{
					for (int i = n - 1; i >= 0; --i)
					{
						T erv;
						int d = FastTest(m, x, mx, v, i, step, erv);
						if (erv < ErV)
						{
							x[i] += step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}

						d = FastTest(m, x, mx, v, i, -step, erv);
						if (erv < ErV)
						{
							x[i] -= step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}
					}
				}

				if (nDirection == 1)
				{
					for (int i = 0; i < n; ++i)
					{
						T erv;
						int d = FastTest(m, x, mx, v, i, step, erv);
						if (erv < ErV)
						{
							x[i] += step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}

						d = FastTest(m, x, mx, v, i, -step, erv);
						if (erv < ErV)
						{
							x[i] -= step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}
					}
				}

				if (cc == 0)
				{
					++stepIndex;
					if (lastIsOK)
						Quicksort(steps, stepCount, xmv::AdaBoostStepValueCompare<T>);
				}
				else
				{
					steps[stepIndex].usedCount += 1;
					stepIndex = 0;
				}
				if (stepIndex >= stepCount) break;

				step = steps[stepIndex].step;

				lastIsOK = (cc != 0);
				//xTemp 
				////cout << x << endl;
			}
			this->parameters = x;
			////cout << '\n' << x << endl;

			return x;
		}

		Vector<T> GradientdescentLearn(Matrix<T> m, Vector<T> initX, Vector<T>v,
			T stepStart = 0.0005, T stepIncreasement = 0.0005, T stepMax = 1, int maxN = -1, int nDirection = 1)
		{
			T step = stepStart;

			Vector<T> x = initX;
			//SolveApproximateLinearEquations(m, v);

			int n = maxN < 0 ? x.Size() : maxN;
			T ErV = -1;
			for (;;)
			{
				Vector<T> mx;
				T OutErV;
				int Er = Test(m, x, mx, v, OutErV);
				if (ErV == -1 || OutErV <= ErV)
					ErV = OutErV;
				cout << Er << " / " << m.Row() << " | ERV:" << ErV << ", step: " << step << "\r";

				int cc = 0;

				if (nDirection == -1)
				{
					for (int i = n - 1; i >= 0; --i)
					{
						T erv;
						int d = FastTest(m, x, mx, v, i, step, erv);
						if (erv < ErV)
						{
							x[i] += step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}

						d = FastTest(m, x, mx, v, i, -step, erv);
						if (erv < ErV)
						{
							x[i] -= step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}
					}
				}

				if (nDirection == 1)
				{
					for (int i = 0; i < n; ++i)
					{
						T erv;
						int d = FastTest(m, x, mx, v, i, step, erv);
						if (erv < ErV)
						{
							x[i] += step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}

						d = FastTest(m, x, mx, v, i, -step, erv);
						if (erv < ErV)
						{
							x[i] -= step;
							Er = d;
							ErV = erv;
							++cc;
							continue;
						}
					}
				}

				if (cc == 0) step += stepIncreasement;
				else step = stepStart;

				if (step > stepMax) break;
				//xTemp 
				////cout << x << endl;
			}
			this->parameters = x;
			////cout << '\n' << x << endl;

			return x;
		}



	public:
		//AdaBoost学习
		//filename，学习后数据保存路径
		//stepStart：梯度下降步长
		//stepIncreasement：题目下降递增
		//stepMax：最大梯度
		//errIncreasement: 对于错误样本的惩罚系数
		//maxCount：最大学习次数
		void AdaBoostLearning(string filename, T stepStart = 0.0005, T stepIncreasement = 0.0005, T stepMax = 0.5,
			T errIncreasement = 1, int maxCount = 1000)
		{
			int stepCount = (stepMax - stepStart) / stepIncreasement + 1;
			AdaBoostStepValue<T> * tarr = NEW("ADABL:TARR") AdaBoostStepValue<T>[stepCount];

			for (int i = 0; i < stepCount; ++i)
			{
				tarr[i].step = stepStart + stepIncreasement * i;
				tarr[i].usedCount = 0;
			}
			Matrix<T> m;
			Vector<T> v;
			T errValue;
			GetSampleAndResultMatrix(m, v);
			Vector<T> v0 = v;

			//cout << m.Col() << endl;

			abParameters = Matrix<T>(m.Col(), maxCount, 0);
			abPowers = Vector<T>(maxCount, 0);

			for (int id = 0; id < maxCount; ++id)
			{
				cout << format(3, 0) << "-------------------------------" << endl;
				Vector<T> initX = SolveApproximateLinearEquations(m, v);
				Vector<T> x = GradientdescentLearn(m, initX, v, tarr, stepCount);
				Quicksort(tarr, stepCount, xmv::AdaBoostStepValueCompare<T>);
				Vector<T> mx;
				int Err = Test(m, x, mx, v, errValue, errIncreasement);

				cout << id << " " << Err << "/" << m.Row() << " ";

				abParameters.SetCol(id, x);
				double k = 1.0 - Err * 1.0 / m.Row();
				cout << k * 100 << "% ";
				if (k > 0.99999) k = 0.99999;
				if (k < 0.00001) k = 0.00001;
				abPowers[id] = std::log(k / (1 - k));

				//检查
				Matrix<T> r = m * abParameters;

				T * ptr;
				enumArray(r.Buffer(), 0, r.Size(), ptr)
					*ptr = ((*ptr >= 0) ? 1 : -1);
				enumArrayEnd;
				Vector<T> RP;
				int errCount = Test(r, abPowers, RP, v0, errValue, 0);

				cout << errCount << " " << (1.0 - errCount * 1.0 / m.Row()) * 100 << "%               \n";

				abPowers = GradientdescentLearn(r, abPowers, v0, 0.003, 0.003, 0.25, id + 1, -1);
				//abPowers.SetMinTo(0);

				errCount = Test(r, abPowers, RP, v0, errValue, 0);

				cout << "                      " << errCount << " " << (1.0 - errCount * 1.0 / m.Row()) * 100 << "%";
				cout << "              " << endl;

				SaveAdaBoostParameters(filename);
			}

			DELETEARR tarr;
		}

		void SaveAdaBoostParameters(string filename)
		{
			ofstream o(filename);
			o << abParameters.Row() << endl;
			o << abParameters.Col() << endl;
			o << format(8, 10) << abParameters << endl;
			o << abPowers.Size() << endl;
			o << format(8, 10) << abPowers << endl;
			o.close();
		}

		void LoadAdaBoostParameters(string filename)
		{
			ifstream o(filename);
			int abpR, abpC;
			o >> abpR >> abpC;
			abParameters = Matrix<T>(abpR, abpC);
			o >> abParameters;

			int ps;
			o >> ps;
			abPowers = Vector<T>(ps);
			o >> abPowers;

			o.close();

			//去掉0
			int end = abPowers.Size();
			for (int i = abPowers.Size() - 1; i >= 0; --i)
			{
				if (Abs<T>(abPowers[i]) > epsilon)
				{
					end = i + 1;
					break;
				}
			}

			////cout << "end:" << end << endl;

			if (end < abPowers.Size())
			{
				abPowers = abPowers.Cut(0, end);
				abParameters = abParameters.Cut(0, abParameters.Row(), 0, end);
			}


		}

		T AdaBoostExecute(Vector<T> & sample)
		{
			if (this->parameterCount <= 0)
			{
				parameterCount = sample.Size();
				parameters = Vector<T>(parameterCount * power + (mul ? parameterCount : 0) + 1, (T)0);
			}

			////cout << "ok" << endl;
			int len = this->ParameterSize();
			////cout << len << endl;
			Vector<T> s(len);
			s[len - 1] = 1;
			T * sp, *sa;
			T first = sample[0];
			T * lineEnd = sample.Buffer() + sample.Size();
			enumArray2(sample.Buffer(), s.Buffer(), 0, sample.Size(), sa, sp)
			{
				T v = (T)1;
				for (int k = 1; k <= power; ++k)
					*sp++ = v *= *sa;

				if (mul)
				{
					if (sa + 1 < lineEnd)
						*sp++ = *sa * sa[1];
					else
						*sp++ = *sa * first;
				}

				--sp;
			}
			enumArrayEnd;

			////cout << s.Size() << " * " << abParameters.Row() << "*" << abParameters.Col() << endl;
			VectorH<T> r = s.Transpose() * abParameters;
			T * ptr;

			enumArray(r.Buffer(), 0, r.Size(), ptr)
			{
				*ptr -= move; // ((*ptr >= move) ? 1 : -1);
				if (*ptr < -1) *ptr = -1;
				if (*ptr > 1) *ptr = 1;
			}
			enumArrayEnd;

			////cout << r.Size() << " * " << abPowers.Size() << endl;
			////cout << r << endl;
			return r * abPowers;
			////cout << "Ada Boost Over ...\r";
		}

		T Execute(Vector<T> & sample)
		{
			if (this->parameterCount <= 0)
			{
				parameterCount = sample.Size();
				parameters = Vector<T>(parameterCount * power + (mul ? parameterCount : 0) + 1, (T)0);
			}

			int len = this->ParameterSize();
			Vector<T> s(len);
			s[len - 1] = 1;
			T * sp, *sa;
			T first = sample[0];
			T * lineEnd = sample.Buffer() + sample.Size();
			enumArray2(sample.Buffer(), s.Buffer(), 0, sample.Size(), sa, sp)
			{
				T v = (T)1;
				for (int k = 1; k <= power; ++k)
					*sp++ = v *= *sa;

				if (mul)
				{
					if (sa + 1 < lineEnd)
						*sp++ = *sa * sa[1];
					else
						*sp++ = *sa * first;
				}

				--sp;
			}
			enumArrayEnd;

			return s * parameters;

		}

		void SaveParameters(string path)
		{
			ofstream o(path);
			o << this->parameters;
			o.close();
		}

		void LoadParameters(string path)
		{
			ifstream i(path);
			i >> this->parameters;
			i.close();
		}
	};
}