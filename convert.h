
#pragma once

#include <common.h>
#include <string>
#include <cstring>
#include <list>
#include <fstream>
#include <iostream>
#include <iomanip>
#ifndef NOSTLFS
#include <filesystem>
#endif
#include <array>

using namespace std;

typedef unsigned char byte;

namespace xmv
{

	string operator + (string & str, int value);
	string operator + (string & str, float value);
	string operator + (string & str, double value);
	string operator + (int value   , string & str);
	string operator + (float value , string & str);
	string operator + (double value, string & str);
	string & operator += (string & str, int value);
	string & operator += (string & str, float value);
	string & operator += (string & str, double value);


#ifdef CLR
	string ToString(System::String ^ str, System::Text::Encoding ^ encoding = System::Text::Encoding::ASCII);
	System::String ^ ToString(string str, System::Text::Encoding ^ encoding = System::Text::Encoding::ASCII);
#endif

	template<class ListOfString>
	void Split(string str, ListOfString & list, char splitChar)
	{
		int len = str.size();
		char * strArr = new char[len + 1];
		strcpy_s(strArr, len + 1, str.c_str());
		char * strArrEnd = strArr + len;
		for(char * p = strArr; p < strArrEnd; ++p)
			if (*p == splitChar) *p = 0;

		for(char * p = strArr; p < strArrEnd; )
		{
			while(!(*p) && p < strArrEnd) ++p;
			if (p >= strArrEnd) break;

			if (*p) list.push_back(p);
			while(*p++);
		}
	}
}