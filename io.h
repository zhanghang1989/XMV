#pragma once

#include <string>
#include <cstring>
#include <list>
#include <fstream>
#include <iostream>
#include <iomanip>


#ifndef NOSTLFS
#include <filesystem>
#endif

#include "common.h"

using namespace std;
typedef unsigned char byte;


namespace xmv
{
	class DigtalFormat
	{
	public:
		int space;
		int decimalPlaces;

	public:
		DigtalFormat() : space(8), decimalPlaces(4) { }
	};

	DigtalFormat format(int decimalPlaces, int space);
	ostream & operator << (ostream & o, DigtalFormat & format);

	void WriteString(ostream & o, string val);
	void WriteString(ostream & o, char * val);
	void WriteChar(ostream & o, char val);
	void WriteByte(ostream & o, byte val);
	void WriteInt8(ostream & o, __int8 val);
	void WriteInt16(ostream & o, __int16 val);
	void WriteInt32(ostream & o, __int32 val);
	void WriteInt64(ostream & o, __int64 val);
	void WriteUInt8(ostream & o, unsigned __int8 val);
	void WriteUInt16(ostream & o, unsigned __int16 val);
	void WriteUInt32(ostream & o, unsigned __int32 val);
	void WriteUInt64(ostream & o, unsigned __int64 val);
	void WriteBytes(ostream & o, byte * buffer, int size);
	void WriteZero(ostream & o, int size);

	string ReadString(istream & i, int size);
	void ReadString(istream & i, char * str, int size);
	char ReadChar(istream & i);
	byte ReadByte(istream & i);
	__int8 ReadInt8(istream & i);
	__int16 ReadInt16(istream & i, bool highByteFirst = false);
	__int32 ReadInt32(istream & i, bool highByteFirst = false);
	__int64 ReadInt64(istream & i, bool highByteFirst = false);
	unsigned __int8 ReadUInt8(istream & i);
	unsigned __int16 ReadUInt16(istream & i, bool highByteFirst = false);
	unsigned __int32 ReadUInt32(istream & i, bool highByteFirst = false);
	unsigned __int64 ReadUInt64(istream & i, bool highByteFirst = false);
	void ReadBytes(istream & i, byte * buffer, int size);
	void ReadZero(istream & i, int size);

	float ReadFloat(istream & i, bool highByteFirst);
	double ReadDouble(istream & i, bool highByteFirst);
}