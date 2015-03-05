#include "convert.h"
#include "memory_.h"
#include "common.h"
#include "time_.h"
#include "io.h"
#include "filesys.h"
#include <cstdlib>
#include <string>
#ifndef NOSTLFS
#include <filesystem>
using namespace std;
using namespace std::tr2::sys;
#endif
#include <vector>

namespace xmv
{
	string operator + (string & str, int value)
	{
		char sValue[20];
		_itoa_s<20>(value, sValue, 10);
		return str + sValue;
	}

	string operator + (string & str, float value)
	{
		char sValue[20];
		sprintf_s<20>(sValue, "%f", (double)value);
		return str + sValue;
	}

	string operator + (string & str, double value)
	{
		char sValue[20];
		sprintf_s<20>(sValue, "%f", value);
		return str + sValue;
	}

	string operator + (int value, string & str)
	{
		char sValue[20];
		_itoa_s<20>(value, sValue, 10);
		return sValue + str;
	}

	string operator + (float value, string & str)
	{
		char sValue[20];
		sprintf_s<20>(sValue, "%f", (double)value);
		return sValue + str;
	}

	string operator + (double value, string & str)
	{
		char sValue[20];
		sprintf_s<20>(sValue, "%f", value);
		return sValue + str;
	}

	string & operator += (string & str, int value)
	{
		char sValue[20];
		_itoa_s<20>(value, sValue, 10);
		return str += sValue;
	}

	string & operator += (string & str, float value)
	{
		char sValue[20];
		sprintf_s<20>(sValue, "%f", (double)value);
		return str += sValue;
	}

	string & operator += (string & str, double value)
	{
		char sValue[20];
		sprintf_s<20>(sValue, "%f", value);
		return str += sValue;
	}

	MemoryChecker MemoryChecker::memoryChecker;
	bool ____randomizeOk = false;
	DigtalFormat ____currentDigtalFormat;
	double epsilon = EPS;

	DigtalFormat format(int decimalPlaces, int space)
	{
		DigtalFormat f;
		f.decimalPlaces = decimalPlaces;
		f.space = space;
		return f;
	}

	ostream & operator << (ostream & o, DigtalFormat & format)
	{
		____currentDigtalFormat = format;
		o << fixed << setprecision(____currentDigtalFormat.decimalPlaces);
		return o;
	}

	void WriteString(ostream & o, string val)
	{
		o << val;
	}

	void WriteString(ostream & o, char * val)
	{
		o << val;
	}

	void WriteChar(ostream & o, char val)
	{
		o << val;
	}

	void WriteByte(ostream & o, byte val)
	{
		o.write((char *)&val, 1);
	}

	void WriteInt8(ostream & o, __int8 val)
	{
		o.write((char *)&val, 1);
	}

	void WriteInt16(ostream & o, __int16 val)
	{
		o.write((char *)&val, 2);
	}

	void WriteInt32(ostream & o, __int32 val)
	{
		o.write((char *)&val, 4);
	}

	void WriteInt64(ostream & o, __int64 val)
	{
		o.write((char *)&val, 8);
	}

	void WriteUInt8(ostream & o, unsigned __int8 val)
	{
		o.write((char *)&val, 1);
	}

	void WriteUInt16(ostream & o, unsigned __int16 val)
	{
		o.write((char *)&val, 2);
	}

	void WriteUInt32(ostream & o, unsigned __int32 val)
	{
		o.write((char *)&val, 4);
	}

	void WriteUInt64(ostream & o, unsigned __int64 val)
	{
		o.write((char *)&val, 8);
	}

	void WriteBytes(ostream & o, byte * buffer, int size)
	{
		o.write((char *)buffer, size);
	}

	void WriteZero(ostream & o, int size)
	{
		char c = 0;
		for (int i = 0; i < size; ++i)
			o.write(&c, 1);
	}

	string ReadString(istream & i, int size)
	{
		char * str = new char[size + 1];
		i.read(str, size);
		str[size] = 0;
		string s(str);
		delete[]str;
		return s;
	}

	void ReadString(istream & i, char * str, int size)
	{
		i.read(str, size);
	}

	char ReadChar(istream & i)
	{
		char c;
		i.read((char *)(&c), 1);
		return c;
	}

	byte ReadByte(istream & i)
	{
		byte c;
		i.read((char *)(&c), 1);
		return c;
	}

	__int8 ReadInt8(istream & i)
	{
		__int8 c;
		i.read((char *)(&c), 1);
		return c;
	}

	void ReadAndFillBufferHighByteFirst(istream & i, void * buffer, int size)
	{
		char * buf = (char *)buffer;
		buf += size - 1;
		for (int j = size - 1; j >= 0; --j, --buf)
		{
			i.read(buf, 1);
		}
	}

	__int16 ReadInt16(istream & i, bool highByteFirst)
	{
		__int16 c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 2);
		else
			i.read((char *)(&c), 2);
		return c;
	}

	__int32 ReadInt32(istream & i, bool highByteFirst)
	{
		__int32 c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 4);
		else
			i.read((char *)(&c), 4);
		return c;
	}

	__int64 ReadInt64(istream & i, bool highByteFirst)
	{
		__int64 c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 8);
		else
			i.read((char *)(&c), 8);
		return c;
	}

	unsigned __int8 ReadUInt8(istream & i)
	{
		unsigned __int8 c;
		i.read((char *)(&c), 1);
		return c;
	}

	unsigned __int16 ReadUInt16(istream & i, bool highByteFirst)
	{
		unsigned __int16 c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 2);
		else
			i.read((char *)(&c), 2);
		return c;
	}

	unsigned __int32 ReadUInt32(istream & i, bool highByteFirst)
	{
		unsigned __int32 c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 4);
		else
			i.read((char *)(&c), 4);
		return c;
	}

	unsigned __int64 ReadUInt64(istream & i, bool highByteFirst)
	{
		unsigned __int64 c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 8);
		else
			i.read((char *)(&c), 8);
		return c;
	}

	void ReadBytes(istream & i, byte * buffer, int size)
	{
		i.read((char *)buffer, size);
	}

	void ReadZero(istream & i, int size)
	{
		char * str = new char[size];
		i.read(str, size);
		delete[]str;
	}

	float ReadFloat(istream & i, bool highByteFirst)
	{
		float c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 4);
		else
			i.read((char *)(&c), 4);
		return c;
	}

	double ReadDouble(istream & i, bool highByteFirst)
	{
		double c;
		if (highByteFirst)
			ReadAndFillBufferHighByteFirst(i, &c, 8);
		else
			i.read((char *)(&c), 4);
		return c;
	}

#ifdef CLR
	string ToString(System::String ^ str, System::Text::Encoding ^ encoding)
	{
		cli::array<byte, 1> ^ bytes = encoding->GetBytes(str);
		char * s = new char[bytes->Length + 1];
		for (int i = 0; i < bytes->Length; ++i)
		{
			s[i] = bytes[i];
		}
		s[bytes->Length] = 0;
		string st = s;
		delete[]s;
		return st;
	}

	System::String ^ ToString(string str, System::Text::Encoding ^ encoding)
	{
		System::String ^ string = gcnew System::String(str.c_str());
		return string;
	}
#endif

#ifndef NOSTLFS
	void CreateDirectory(string directoryName)
	{
		path dname(directoryName);
		create_directories(dname);
	}

	string StartupPath()
	{
		return initial_path<path>();
	}
#endif


	ostream & operator << (ostream & o, xmv::DateTime dt)
	{
		o << dec;
		o << dt.year << "-";
		if (dt.month < 10) o << "0";
		o << dt.month << "-";
		if (dt.day < 10) o << "0";
		o << dt.day << " ";

		if (dt.hour < 10) o << "0";
		o << dt.hour << ":";
		if (dt.minute < 10) o << "0";
		o << dt.minute << ":";
		if (dt.second < 10) o << "0";
		o << dt.second << ".";

		if (dt.millisecond < 100) o << "0";
		if (dt.millisecond < 10) o << "0";
		o << (int)dt.millisecond;

		return o;
	}

	ostream & operator << (ostream & o, xmv::TimeSpan dt)
	{
		if (dt.day < 10) o << "0";
		o << dt.day << ".";

		if (dt.hour < 10) o << "0";
		o << dt.hour << ":";
		if (dt.minute < 10) o << "0";
		o << dt.minute << ":";
		if (dt.second < 10) o << "0";
		o << dt.second << ".";

		if (dt.millisecond < 100) o << "0";
		if (dt.millisecond < 10) o << "0";
		o << (int)dt.millisecond;

		return o;
	}

	DateTime operator + (DateTime dt, TimeSpan span)
	{
		return DateTime(dt.ticks + span.ticks);
	}

	DateTime operator - (DateTime dt, TimeSpan span)
	{
		return DateTime(dt.ticks - span.ticks);
	}

	TimeSpan operator - (DateTime dt1, DateTime dt2)
	{
		return TimeSpan(dt1.ticks - dt2.ticks);
	}

	TimeSpan operator + (TimeSpan dt1, TimeSpan dt2)
	{
		return TimeSpan(dt1.ticks + dt2.ticks);
	}

	bool operator > (TimeSpan dt1, TimeSpan dt2)
	{
		return dt1.ticks > dt2.ticks;
	}

	bool operator < (TimeSpan dt1, TimeSpan dt2)
	{
		return dt1.ticks < dt2.ticks;
	}

	bool operator >= (TimeSpan dt1, TimeSpan dt2)
	{
		return dt1.ticks >= dt2.ticks;
	}

	bool operator <= (TimeSpan dt1, TimeSpan dt2)
	{
		return dt1.ticks <= dt2.ticks;
	}

	bool operator == (TimeSpan dt1, TimeSpan dt2)
	{
		return dt1.ticks == dt2.ticks;
	}

	bool operator != (TimeSpan dt1, TimeSpan dt2)
	{
		return dt1.ticks != dt2.ticks;
	}

	bool operator > (DateTime dt1, DateTime dt2)
	{
		return dt1.ticks > dt2.ticks;
	}

	bool operator < (DateTime dt1, DateTime dt2)
	{
		return dt1.ticks < dt2.ticks;
	}

	bool operator >= (DateTime dt1, DateTime dt2)
	{
		return dt1.ticks >= dt2.ticks;
	}

	bool operator <= (DateTime dt1, DateTime dt2)
	{
		return dt1.ticks <= dt2.ticks;
	}

	bool operator == (DateTime dt1, DateTime dt2)
	{
		return dt1.ticks == dt2.ticks;
	}

	bool operator != (DateTime dt1, DateTime dt2)
	{
		return dt1.ticks != dt2.ticks;
	}

	string toLower(const string & str)
	{
		char * buff = new char[str.length() + 1];
		char * p2 = buff;
		for (const char * p = str.c_str(); *p; ++p, ++p2)
		{
			if (*p >= 'A' && *p <= 'Z')
				*p2 = *p + 32;
			else
				*p2 = *p;
		}
		*p2 = 0;
		string s(buff);
		delete buff;
		return s;
	}

	string GetFileName(string path)
	{
		int k = path.find_last_of('\\');
		if (k >= 0)
			return path.substr(k + 1);
		else
			return path;
	}
}