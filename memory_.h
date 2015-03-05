#pragma once
#include <cstdlib>
#include <list>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

namespace xmv
{
	struct MemoryItem
	{
		string description;
		void * address;
	};

	class MemoryChecker
	{
	private:
		list<MemoryItem *> * items;

		string tempDescription;
		MemoryChecker() : items(NULL)
		{
		}

		~MemoryChecker()
		{
			if (items) delete items;
			items = NULL;
		}

	public:
		static MemoryChecker memoryChecker;

		template<typename T> T* operator + (T* address)
		{
			if (items == NULL) items = new list<MemoryItem *>;
			MemoryItem * item = new MemoryItem;
			item->description = tempDescription;
			item->address = address;
			items->push_front(item);
			return address;
		}

		MemoryChecker & operator += (string description)
		{
			this->tempDescription = description;
			return *this;
		}

		template<typename T> MemoryChecker & operator - (T * address)
		{
			if (! items)
			{
				delete address;
				return *this;
			}

			list<MemoryItem *>::iterator end = items->end();
			for(list<MemoryItem *>::iterator i = items->begin(); i != end; ++i)
			{
				MemoryItem * ip = *i;
				if (ip->address == address)
				{
					delete address;
					delete ip;
					items->erase(i);
					return *this;
				}
			}
			throw;
		}

		template<typename T> MemoryChecker & operator -= (T * address)
		{
			if (! items)
			{
				delete []address;
				return *this;
			}
			list<MemoryItem *>::iterator end = items->end();
			for(list<MemoryItem *>::iterator i = items->begin(); i != end; ++i)
			{
				MemoryItem * ip = *i;
				if (ip->address == address)
				{
					if (!address) return *this; //
					delete[]address;
					if (ip) delete ip;
					items->erase(i);
					return *this;
				}
			}
			delete []address;
			
			return *this;
		}

		void Output(string filename)
		{
#if (defined DEBUG) || (defined _DEBUG)
			ofstream o(filename);
			list<MemoryItem*>::iterator end = items->end();
			for(list<MemoryItem*>::iterator i = items->begin(); i != end; ++i)
			{
				MemoryItem * ip = *i;
				o << ip->address << ": " << ip->description << "\r\n";
			}
#endif
		}
	};


#if (defined DEBUG) || (defined _DEBUG)
#define NEW(description) (xmv::MemoryChecker::memoryChecker += (description)) + new 
#define DELETE xmv::MemoryChecker::memoryChecker - 
#define DELETEARR xmv::MemoryChecker::memoryChecker -= 
#else
#define NEW(description) new 
#define DELETE delete 
#define DELETEARR delete [] 
#endif


}