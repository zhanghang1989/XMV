#pragma once

#include "common.h"
#include <cstdlib>
#include <iostream>
#include <cstdio>
using namespace std;

namespace xmv
{
	const unsigned long long ticksOfMicroseconds = 10;
	const unsigned long long ticksOfMillisecond = ticksOfMicroseconds * 1000;
	const unsigned long long ticksOfSecond = ticksOfMillisecond * 1000;
	const unsigned long long ticksOfMinute = ticksOfSecond * 60;
	const unsigned long long ticksOfHour = ticksOfMinute * 60;
	const unsigned long long ticksOfDay = ticksOfHour * 24;
	const int daysOfYear = 365;
	const int daysOfLeapYear = 366;
	const int daysOf400Years = daysOfYear * 400 + 97;

	class TimeSpan;
	class DateTime
	{
	private: 
		unsigned long long ticks;
		int year;
		int month;
		int day;
		int hour;
		int minute;
		int second;
		double millisecond;

	public: unsigned long long GetTicks() { return ticks; }
	public: int GetYear() { return year;  }
	public: int GetMonth() { return month;  }
	public: int GetDay() { return day;  }
	public: int GetHour() { return hour;  }
	public: int GetMinute() { return minute;  }
	public: int GetSecond() { return second;  }
	public: double Millisecond() { return millisecond;  }

	public : static bool IsLeapYear(int year)
			 {
				 if (year % 400 == 0) return true;
				 if (year % 100 == 0) return false;
				 if (year % 4 == 0) return true;
				 return false;
			 }

	public: DateTime()
			{
				this->year = 1;
				this->month = 1;
				this->day = 1;
				this->hour = 0;
				this->minute = 0;
				this->second = 0;
				this->millisecond = 0;
				this->ticks = 0;
			}

	public: DateTime(int year, int month, int day, int hour = 0, int minute = 0, int second = 0, double millisecond = 0)
			{
				this->year = year;
				this->month = month;
				this->day = day;
				this->hour = hour;
				this->minute = minute;
				this->second = second;
				this->millisecond = millisecond;

				int days = year / 400 * daysOf400Years;
				int y = year / 400 * 400 + 1;
				for(; y < year; ++y)
				{
					days += IsLeapYear(y) ? daysOfLeapYear : daysOfYear;
				}

				int daysInMonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
				if (IsLeapYear(year)) daysInMonth[1] = 29;
				for(int m = 1; m < month; ++m)
					days += daysInMonth[m - 1];
				days += day - 1;

				ticks = days * ticksOfDay +
					hour * ticksOfHour +
					minute * ticksOfMinute +
					second * ticksOfSecond +
					(unsigned long long)(millisecond * ticksOfMillisecond);
			}

	public: DateTime(unsigned long long ticks)
			{
				this->ticks = ticks;

				unsigned long long days = ticks / ticksOfDay;

				year = days / daysOf400Years * 400 + 1;
				days %= daysOf400Years;
				for( ; ; year++)
				{
					int dayInThisYear = IsLeapYear(year) ? daysOfLeapYear : daysOfYear;
					if (days > dayInThisYear) 
					{
						days -= dayInThisYear;
					}
					else break;
				}

				int daysInMonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
				if (IsLeapYear(year)) daysInMonth[1] = 29;
				for(month = 1; month <= 12; month++)
				{
					if (days >= daysInMonth[month - 1])
					{
						days -= daysInMonth[month - 1];
						continue;
					}
					break;
				}

				day = days + 1;

				unsigned long long ticksInDay = ticks % ticksOfDay;
				this->hour = ticksInDay / ticksOfHour;
				this->minute = (ticksInDay %= ticksOfHour) / ticksOfMinute;
				this->second = (ticksInDay %= ticksOfMinute) / ticksOfSecond;
				this->millisecond = (ticksInDay %= ticksOfSecond) / (double)ticksOfMillisecond;
			}

			friend ostream & operator << (ostream & o, DateTime dt);
			friend DateTime operator + (DateTime dt, TimeSpan span);
			friend DateTime operator - (DateTime dt, TimeSpan span);
			friend TimeSpan operator - (DateTime dt1, DateTime dt2);

			friend bool operator > (DateTime dt1, DateTime dt2);
			friend bool operator < (DateTime dt1, DateTime dt2);
			friend bool operator >= (DateTime dt1, DateTime dt2);
			friend bool operator <= (DateTime dt1, DateTime dt2);
			friend bool operator == (DateTime dt1, DateTime dt2);
			friend bool operator != (DateTime dt1, DateTime dt2);
	};


	class TimeSpan
	{
	private: 
		unsigned long long ticks;
		int day;
		int hour;
		int minute;
		int second;
		double millisecond;

	public: unsigned long long GetTicks() { return ticks; }
	public: int GetDay() { return day;  }
	public: int GetHour() { return hour;  }
	public: int GetMinute() { return minute;  }
	public: int GetSecond() { return second;  }
	public: double Millisecond() { return millisecond;  }

	public: TimeSpan()
			{
				this->day = 0;
				this->hour = 0;
				this->minute = 0;
				this->second = 0;
				this->millisecond = 0;
			}

	public: TimeSpan(int day, int hour, int minute , int second, double millisecond)
			{
				this->day = day;
				this->hour = hour;
				this->minute = minute;
				this->second = second;
				this->millisecond = millisecond;

				ticks = day * ticksOfDay +
					hour * ticksOfHour +
					minute * ticksOfMinute +
					second * ticksOfSecond +
					(unsigned long long)(millisecond * ticksOfMillisecond);
			}

	public: TimeSpan(unsigned long long ticks)
			{
				this->ticks = ticks;

				this->day = ticks / ticksOfDay;

				unsigned long long ticksInDay = ticks % ticksOfDay;
				this->hour = ticksInDay / ticksOfHour;
				this->minute = (ticksInDay %= ticksOfHour) / ticksOfMinute;
				this->second = (ticksInDay %= ticksOfMinute) / ticksOfSecond;
				this->millisecond = (ticksInDay %= ticksOfSecond) / (double)ticksOfMillisecond;
			}

			friend ostream & operator << (ostream & o, TimeSpan dt);
			friend DateTime operator + (DateTime dt, TimeSpan span);
			friend DateTime operator - (DateTime dt, TimeSpan span);
			friend TimeSpan operator - (DateTime dt1, DateTime dt2);

			friend TimeSpan operator + (TimeSpan dt1, TimeSpan dt2);

			friend bool operator > (TimeSpan dt1, TimeSpan dt2);
			friend bool operator < (TimeSpan dt1, TimeSpan dt2);
			friend bool operator >= (TimeSpan dt1, TimeSpan dt2);
			friend bool operator <= (TimeSpan dt1, TimeSpan dt2);
			friend bool operator == (TimeSpan dt1, TimeSpan dt2);
			friend bool operator != (TimeSpan dt1, TimeSpan dt2);


	};


	ostream & operator << (ostream & o, DateTime dt);
	ostream & operator << (ostream & o, TimeSpan dt);
	DateTime operator + (DateTime dt, TimeSpan span);
	DateTime operator - (DateTime dt, TimeSpan span);
	TimeSpan operator - (DateTime dt1, DateTime dt2);
	TimeSpan operator + (TimeSpan dt1, TimeSpan dt2);

	bool operator > (TimeSpan dt1, TimeSpan dt2);
	bool operator < (TimeSpan dt1, TimeSpan dt2);
	bool operator >= (TimeSpan dt1, TimeSpan dt2);
	bool operator <= (TimeSpan dt1, TimeSpan dt2);
	bool operator == (TimeSpan dt1, TimeSpan dt2);
	bool operator != (TimeSpan dt1, TimeSpan dt2);

	bool operator > (TimeSpan dt1, TimeSpan dt2);
	bool operator < (TimeSpan dt1, TimeSpan dt2);
	bool operator >= (TimeSpan dt1, TimeSpan dt2);
	bool operator <= (TimeSpan dt1, TimeSpan dt2);
	bool operator == (TimeSpan dt1, TimeSpan dt2);
	bool operator != (TimeSpan dt1, TimeSpan dt2);

}