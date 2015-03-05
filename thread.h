#include <list>

#ifdef CLR

using namespace System::Threading;
using namespace System::Collections::Generic;

namespace xmv
{

	template <class argumentType>
	ref class StartThreadGroupArg
	{
	public :
		argumentType * arg;
		void (* function)(argumentType * arguments);
	};

	template <class argumentType>
	void ThreadRunOnce(System::Object ^ argObj)
	{
		StartThreadGroupArg<argumentType> ^ arg = (StartThreadGroupArg<argumentType> ^)argObj;
		arg->function(arg->arg);
	}

	template <class argumentType>
	void StartThreadGroup(void (* threadFunction)(argumentType * arguments), list<argumentType *> & argList)
	{
		int k1 = System::DateTime::Now.Ticks;
		List<Thread ^> ^ threadList = gcnew List<Thread ^>();
		auto end = argList.end();
		for(auto i = argList.begin(); i != end; ++i)
		{
			ParameterizedThreadStart ^ threadHandler = gcnew ParameterizedThreadStart(ThreadRunOnce<argumentType>);
			StartThreadGroupArg<argumentType> ^ threadArg = gcnew StartThreadGroupArg<argumentType>();
			threadArg->arg = *i;
			threadArg->function = threadFunction;

			Thread ^ thread = gcnew Thread(threadHandler);
			thread->Priority = ThreadPriority::Highest;
			thread->Start(threadArg);

			threadList->Add(thread);
		}
		int k2 = System::DateTime::Now.Ticks;

		//cout << (k2 - k1) / 10000.0 << "ms  " ;

		bool allcompleted = false;
		do
		{
			Thread::Sleep(10);

			allcompleted = true;
			for each (Thread ^ thread in threadList)
			{
				ThreadState ts = thread->ThreadState;
				if (ts != ThreadState::Aborted && ts != ThreadState::Stopped)
				{
					allcompleted = false;
					break;
				}
			}
		}while(!allcompleted);
	}

	template <class argumentType>
	ref class StartThreadRunArg
	{
	public :
		argumentType * arg;
		void (* function)(argumentType * arguments);
		bool starting;
		bool finished;
		bool exiting;
		int id;

		StartThreadRunArg()
		{
			id = 0;
			arg = nullptr;
			function = nullptr;
			starting = false;
			finished = false;
			exiting = false;
		}
	};




	template <class argumentType>
	void ThreadRun(System::Object ^ argObj)
	{
		StartThreadRunArg<argumentType> ^ arg = (StartThreadRunArg<argumentType> ^)argObj;

		while(! arg->exiting)
		{
			if (arg->starting)
			{
				arg->starting = false;
				arg->finished = false;
				arg->function(arg->arg);
				//cout << arg->id << " finished!" << endl;
				arg->finished = true;
			}

			Thread::Sleep(5);
		}
	}

	template <class argumentType>
	ref class ThreadPool
	{
	public: 
		ThreadPool()
		{
			threadList = gcnew List<Thread ^>();
			argList = gcnew List<StartThreadRunArg<argumentType> ^>();
		}

	private:
		List<Thread ^> ^ threadList;
		List<StartThreadRunArg<argumentType> ^> ^ argList;

	public:
		void Run(void (* threadFunction)(argumentType * arguments), int threadCount)
		{
			int k1 = System::DateTime::Now.Ticks;

			for(int i = 0; i < threadCount; ++i)
			{
				ParameterizedThreadStart ^ threadHandler = gcnew ParameterizedThreadStart(ThreadRun<argumentType>);
				StartThreadRunArg<argumentType> ^ threadArg = gcnew StartThreadRunArg<argumentType>();
				threadArg->arg = nullptr;
				threadArg->function = threadFunction;
				threadArg->id = i;

				Thread ^ thread = gcnew Thread(threadHandler);
				thread->Start(threadArg);
				thread->Priority = ThreadPriority::Highest;

				threadList->Add(thread);
				argList->Add(threadArg);
			}
			int k2 = System::DateTime::Now.Ticks;

			cout << (k2 - k1) / 10000.0 << "ms  \n" ;
		}

		void Job(list<argumentType *> & argList)
		{
			bool allcompleted = false;


			auto end = argList.end();
			int id = 0;
			for (auto i = argList.begin(); i != end; ++i, ++id)
			{
				if (id < this->argList->Count)
				{
					StartThreadRunArg<argumentType> ^ arg = this->argList[id];
					arg->arg = *i;
					arg->finished = false;
					arg->starting = true;
				}
			}

			do
			{
				Thread::Sleep(5);
				allcompleted = true;
				for(int id = 0; id < this->argList->Count; ++id)
				{
					if (! this->argList[id]->finished)
					{
						allcompleted = false;
						break;
					}
				}
			}while(!allcompleted);
		}
	};
}
#endif