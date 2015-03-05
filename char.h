#pragma once

#include <cstdlib>
#include <string>

#include "common.h"
#include "convert.h"
#include "pixcel.h"


using namespace std;


namespace xmv
{

	template <typename ImgT>
	void GetCharMartrix_5x7(ImgT * img, char * ch, int chSize, char c)
	{
		switch(c)
		{                   
		case '0': strcpy_s(ch, chSize, "01110 10001 11001 10101 10011 10001 01110"); break;
		case '1': strcpy_s(ch, chSize, "01110 00100 00100 00100 00100 01100 00100"); break;
		case '2': strcpy_s(ch, chSize, "11111 01000 00100 00010 00001 10001 01110"); break;
		case '3': strcpy_s(ch, chSize, "11110 00001 00001 00110 00001 00001 11110"); break;
		case '4': strcpy_s(ch, chSize, "00001 11111 10001 10001 01001 00101 00011"); break;
		case '5': strcpy_s(ch, chSize, "11110 00001 00001 11110 10000 10000 11111"); break;
		case '6': strcpy_s(ch, chSize, "01110 10001 10001 11110 10000 10000 01111"); break;
		case '7': strcpy_s(ch, chSize, "00100 00100 00100 00010 00010 00001 11111"); break;
		case '8': strcpy_s(ch, chSize, "01110 10001 10001 01110 10001 10001 01110"); break;
		case '9': strcpy_s(ch, chSize, "11100 00010 00001 01111 10001 10001 01110"); break;
		case 'A': case 'a': strcpy_s(ch, chSize, "10001 11111 10001 01010 01010 00100 00100"); break;
		case 'B': case 'b': strcpy_s(ch, chSize, "11110 10001 10001 11110 10001 10001 11110"); break;
		case 'C': case 'c': strcpy_s(ch, chSize, "01110 10001 10000 10000 10000 10001 01110"); break;
		case 'D': case 'd': strcpy_s(ch, chSize, "11110 10001 10001 10001 10001 10001 11110"); break;
		case 'E': case 'e': strcpy_s(ch, chSize, "11111 10000 10000 11110 10000 10000 11111"); break;
		case 'F': case 'f': strcpy_s(ch, chSize, "10000 10000 10000 11110 10000 10000 11111"); break;
		case 'G': case 'g': strcpy_s(ch, chSize, "01111 10001 10001 10111 10000 10000 01111"); break;
		case 'H': case 'h': strcpy_s(ch, chSize, "10001 10001 10001 11111 10001 10001 10001"); break;
		case 'I': case 'i': strcpy_s(ch, chSize, "01110 00100 00100 00100 00100 00100 01110"); break;
		case 'J': case 'j': strcpy_s(ch, chSize, "11100 00010 00010 00010 00010 00010 00111"); break;
		case 'K': case 'k': strcpy_s(ch, chSize, "10001 10010 10100 11000 10100 10010 10001"); break;
		case 'L': case 'l': strcpy_s(ch, chSize, "11111 10000 10000 10000 10000 10000 10000"); break;
		case 'M': case 'm': strcpy_s(ch, chSize, "10001 10001 10101 10101 11011 11011 10001"); break;
		case 'N': case 'n': strcpy_s(ch, chSize, "10001 10011 10011 10101 11001 11001 10001"); break;
		case 'O': case 'o': strcpy_s(ch, chSize, "01110 10001 10001 10001 10001 10001 01110"); break;
		case 'P': case 'p': strcpy_s(ch, chSize, "10000 10000 10000 11110 10001 10001 11110"); break;
		case 'Q': case 'q': strcpy_s(ch, chSize, "00011 01110 10101 10001 10001 10001 01110"); break;
		case 'R': case 'r': strcpy_s(ch, chSize, "10001 10010 10100 11110 10001 10001 11110"); break;
		case 'S': case 's': strcpy_s(ch, chSize, "01110 10001 00001 01110 10000 10001 01110"); break;
		case 'T': case 't': strcpy_s(ch, chSize, "00100 00100 00100 00100 00100 00100 11111"); break;
		case 'U': case 'u': strcpy_s(ch, chSize, "01110 10001 10001 10001 10001 10001 10001"); break;
		case 'V': case 'v': strcpy_s(ch, chSize, "00100 00100 01010 01010 10001 10001 10001"); break;
		case 'W': case 'w': strcpy_s(ch, chSize, "10001 11011 11011 10101 10101 10001 10001"); break;
		case 'X': case 'x': strcpy_s(ch, chSize, "10001 10001 01010 00100 01010 10001 10001"); break;
		case 'Y': case 'y': strcpy_s(ch, chSize, "00100 00100 00100 01010 01010 10001 10001"); break;
		case 'Z': case 'z': strcpy_s(ch, chSize, "11111 10000 01000 00100 00010 00001 11111"); break;
		case ':': strcpy_s(ch, chSize, "00000 00100 00000 00000 00000 00100 00000"); break;
		case ';': strcpy_s(ch, chSize, "01000 00100 00000 00000 00000 00100 00000"); break;
		case '\'': strcpy_s(ch, chSize, "00000 00000 00000 00000 00100 00100 00100"); break;
		case '\"': strcpy_s(ch, chSize, "00000 00000 00000 00000 01010 01010 01010"); break;
		case '.': strcpy_s(ch, chSize, "00100 00000 00000 00000 00000 00000 00000"); break;
		case '`': strcpy_s(ch, chSize, "00000 00000 00000 00000 00000 00100 01000"); break;
		case '@': strcpy_s(ch, chSize, "01110 10000 10110 11011 10101 10001 01110"); break;
		case '^': strcpy_s(ch, chSize, "00000 00000 00000 00000 10001 01010 00100"); break;
		case '$': strcpy_s(ch, chSize, "00100 11110 00101 01110 10100 01111 00100"); break;
		case ',': strcpy_s(ch, chSize, "00100 00010 00110 00000 00000 00000 00000"); break;
		case '~': strcpy_s(ch, chSize, "00000 00000 10010 10101 01001 00000 00000"); break;
		case '-': strcpy_s(ch, chSize, "00000 00000 00000 11111 00000 00000 00000"); break;
		case '+': strcpy_s(ch, chSize, "00000 00100 00100 11111 00100 00100 00000"); break;
		case '=': strcpy_s(ch, chSize, "00000 00000 11111 00000 11111 00000 00000"); break;
		case '*': strcpy_s(ch, chSize, "00100 10101 01110 00100 01110 10101 00100"); break;
		case '(': strcpy_s(ch, chSize, "00010 00100 01000 01000 01000 00100 00010"); break;
		case ')': strcpy_s(ch, chSize, "01000 00100 00010 00010 00010 00100 01000"); break;
		case '{': strcpy_s(ch, chSize, "00110 00100 00100 01000 00100 00100 00110"); break;
		case '}': strcpy_s(ch, chSize, "01100 00100 00100 00010 00100 00100 01100"); break;
		case '[': strcpy_s(ch, chSize, "01110 01000 01000 01000 01000 01000 01110"); break;
		case ']': strcpy_s(ch, chSize, "01110 00010 00010 00010 00010 00010 01110"); break;
		case '<': strcpy_s(ch, chSize, "00010 00100 01000 10000 01000 00100 00010"); break;
		case '>': strcpy_s(ch, chSize, "01000 00100 00010 00001 00010 00100 01000"); break;
		case '!': strcpy_s(ch, chSize, "00100 00000 00100 00100 00100 00100 00100"); break;
		case '%': strcpy_s(ch, chSize, "10011 10011 01000 00100 00010 11001 11001"); break;
		case '#': strcpy_s(ch, chSize, "01010 11111 01010 01010 01010 11111 01010"); break;
		case '&': strcpy_s(ch, chSize, "01101 10010 01101 01100 10010 10010 01100"); break;
		case '_': strcpy_s(ch, chSize, "11111 00000 00000 00000 00000 00000 00000"); break;
		case '?': strcpy_s(ch, chSize, "00100 00000 00100 00010 00001 10001 01110"); break;
		case '/': strcpy_s(ch, chSize, "10000 10000 01000 00100 00010 00001 00001"); break;
		case '\\': strcpy_s(ch, chSize, "00001 00001 00010 00100 01000 10000 10000"); break;
		case '|': strcpy_s(ch, chSize, "00100 00100 00100 00100 00100 00100 00100"); break;
		default:  strcpy_s(ch, chSize, "00000 00000 00000 00000 00000 00000 00000"); break;
		}
	}

	template <typename ImgT>
	void GetCharMartrix_16x18(ImgT * img, char * ch, int chSize, char c)
	{
		switch (c)
		{
		case 'A':
		case 'a': 
			strcpy_s(ch, chSize, "1111000000001111 1111000000001111 1111000000001111 0111100000011110 0111111111111110 0111111111111110 0011110000111100 0011110000111100 0011110000111100 0001111001111000 0001111001111000 0001111001111000 0000111111110000 0000111111110000 0000111111110000 0000011111100000 0000011111100000 0000011111100000"); break;
		case 'B':
		case 'b': 
			strcpy_s(ch, chSize, "1111111111110000 1111111111111100 1111000000111110 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000111110 1111111111111100 1111111111111000 1111000001111100 1111000000011110 1111000000011110 1111000000011110 1111000000011110 1111000001111100 1111111111111000 111111111110000"); break;
		case 'C':
		case 'c': 
			strcpy_s(ch, chSize, "0000011111110000 0001111111111100 0011110000011110 0111100000001111 0111100000001111 1111000000000000 1111000000000000 1111000000000000 1111000000000000 1111000000000000 1111000000000000 1111000000000000 1111000000000000 0111100000001111 0111100000001111 0011110000011110 0001111111111100 0000011111110000"); break;
		case 'D':
		case 'd':
			strcpy_s(ch, chSize, "1111111111100000 1111111111111000 1111000001111100 1111000000011110 1111000000011110 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000011110 1111000000011110 1111000001111100 1111111111111000 1111111111100000"); break;
		case 'X':
		case 'x':
			strcpy_s(ch, chSize, "1111000000001111 0111100000011110 0111100000011110 0011110000111100 0001111001111000 0001111001111000 0000111111110000 0000011111100000 0000011111100000 0000011111100000 0000011111100000 0000111111110000 0001111001111000 0001111001111000 0011110000111100 0111100000011110 1111000000001111 1111000000001111"); break;
		case 'Y':
		case 'y':
			strcpy_s(ch, chSize, "0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000011111100000 0000011111100000 0000111111110000 0000111111110000 0001111001111000 0001111001111000 0011110000111100 0011110000111100 0111100000011110 0111100000011110 1111000000001111 1111000000001111"); break;
		case '0':
			strcpy_s(ch, chSize, "0000011111100000 0001111111111000 0011111001111100 0111100000011110 0111100000011110 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 0111100000011110 0111100000011110 0011111001111100 0001111111111000 0000011111100000"); break;
		case '1':
			strcpy_s(ch, chSize, "0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0001111111000000 0001111111000000 0000111111000000 0000011111000000 0000001111000000"); break;
		case '2':
			strcpy_s(ch, chSize, "1111111111111111 1111111111111111 1111100000000000 0111100000000000 0011110000000000 0001111000000000 0000011110000000 0000000111100000 0000000001111000 0000000000111100 0000000000011110 0000000000001111 0000000000001111 1111000000001111 1111000000001111 0111110000111110 0011111111111100 0000111111110000"); break;
		case '3':
			strcpy_s(ch, chSize, "0000111111110000 0011111111111100 0111110000111110 1111000000001111 1111000000001111 0000000000001111 0000000000011110 0000000001111100 0000001111110000 0000001111110000 0000000001111100 0000000000011110 0000000000001111 1111000000001111 1111000000001111 0111110000111110 0011111111111100 0000111111110000"); break;
		case '4':
			strcpy_s(ch, chSize, "0000000001111000 0000000001111000 0000000001111000 0000000001111000 1111111111111111 1111111111111111 0111100001111000 0111100001111000 0011110001111000 0011110001111000 0001111001111000 0001111001111000 0000111101111000 0000111101111000 0000011111111000 0000011111111000 0000001111111000 0000001111111000"); break;
		case '5':
			strcpy_s(ch, chSize, "0000111111110000 0011111111111100 0111110000111110 1111000000001111 1111000000001111 0000000000001111 0000000000001111 0000000000001111 0000000000001111 1111100000111110 1111111111111100 1111111111110000 1111000000000000 1111000000000000 1111000000000000 1111000000000000 1111111111111111 1111111111111111");  break;
		case '6':
			strcpy_s(ch, chSize, "0000111111110000 0011111111111100 0111110000111110 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111100000001111 1111111000111110 1111111111111100 1111011111110000 1111000000000000 1111000000000000 1111000000001111 0111110000111110 0011111111111100 0000111111110000");  break;
		case '7':
			strcpy_s(ch, chSize, "0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000001111000000 0000000111100000 0000000111100000 0000000111100000 0000000011110000 0000000011110000 0000000001111000 0000000001111000 0000000000111100 0000000000011110 1111111111111111 1111111111111111");  break;
		case '8':
			strcpy_s(ch, chSize, "0000111111110000 0011111111111100 0111110000111110 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 0111110000111110 0011111111111100 0011111111111100 0111110000111110 1111000000001111 1111000000001111 1111000000001111 0111110000111110 0011111111111100 0000111111110000");  break;
		case '9':
			strcpy_s(ch, chSize, "0000111111110000 0011111111111100 0111110000111110 1111000000001111 0000000000001111 0000000000001111 0000111111101111 0011111111111111 0111110001111111 1111000000011111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 1111000000001111 0111110000111110 0011111111111100 0000111111110000");  break;
		default: 
			strcpy_s(ch, chSize, "0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000 0000000000000000"); break;
		}
	}

	enum CharSizeType 
	{
		_5x7,
		_16x18,
	};

	template<typename ImgT>
	int GetCharWidth(ImgT *img,  CharSizeType size = CharSizeType::_5x7)
	{
		switch(size)
		{
		case CharSizeType::_5x7: return 5;
		case CharSizeType::_16x18: return 16;
		}
		return 5;
	}


	template<typename ImgT>
	int GetCharHeight(ImgT *img,  CharSizeType size = CharSizeType::_5x7)
	{
		switch(size)
		{
		case CharSizeType::_5x7: return 7;
		case CharSizeType::_16x18: return 18;
		}
		return 5;
	}

	template <typename ImgT, typename TValue>
	void DrawChar(ImgT * img, int width, int height, int x, int y, TValue value, char c, 
		CharSizeType size = CharSizeType::_5x7, int times = 1)
	{
		int cw = GetCharWidth(img, size);
		int ch = GetCharHeight(img, size);
		int chSize = (cw + 1) * ch + 1;
		char * arr = NEW("DRAW CHAR BUFFER") char[chSize];

		switch(size)
		{
		case CharSizeType::_5x7:
			GetCharMartrix_5x7(img, arr, chSize, c);
			break;
		case CharSizeType::_16x18:
			GetCharMartrix_16x18(img, arr, chSize, c);
			break;
		}

		ImgT * p;
		char * a = arr;
		if (times <= 1)
		{
			enumImage(img, width, height, x, x + cw, y, y + ch, p)
			{
				if (*a++ == '1')
					*p = (ImgT)value;
				if (x_____ == x + cw - 1) ++a;
			} enumImageEnd;
		}
		else
		{
			int x_, y_;
			enumImageStep(img, width, height, x, x + cw * times, y, y + ch * times, times, times, p)
			{
				if (*a++ == '1')
				{
					for(y_ = 0; y_ < times; ++y_)
						for(x_ = 0; x_ < times; ++x_)
						{
							p[x_ + y_ * width] = (ImgT)value;
						}
				}
				if (x_____ == x + cw * times - times) ++a;
			} enumImageEnd;
		}


		DELETEARR(arr);
	}

	template <typename ImgT, typename TValue>
	void DrawString(ImgT * img, int width, int height, int x, int y, TValue value, string text, CharSizeType size = CharSizeType::_5x7, int times = 1)
	{
		DrawString(img, width, height, x, y, value, (char * )text.c_str(), size, times);
	}

	template <typename ImgT, typename TValue>
	void DrawString(ImgT * img, int width, int height, int x, int y, TValue value, TValue borderValue, string text, int borderWidth = 1, CharSizeType size = CharSizeType::_5x7, int times = 1)
	{
		DrawString(img, width, height, x, y, value, borderValue, text.c_str(), borderWidth, size, times);
	}

	template <typename ImgT, typename TValue>
	void DrawString(ImgT * img, int width, int height, int x, int y, TValue value, char * text, CharSizeType size = CharSizeType::_5x7, int times = 1)
	{
		int w = GetCharWidth(img, size);

		for(char * p = text; *p; ++p)
		{
			DrawChar(img, width, height, x, y, value, *p, size, times);
			x += (w + w / 5) * times;
		}
	}

	template <typename ImgT, typename TValue>
	void DrawString(ImgT * img, int width, int height, int x, int y, TValue value, int text, CharSizeType size = CharSizeType::_5x7, int times = 1)
	{
		char v[32];
		_itoa(text, v, 10);
		DrawString(img, width, height, x, y, value, v, size, times);
	}

	template <typename ImgT, typename TValue>
	void DrawString(ImgT * img, int width, int height, int x, int y, TValue value, TValue borderValue, int text, int borderWidth = 1, CharSizeType size = CharSizeType::_5x7, int times = 1)
	{
		char v[32];
		_itoa(text, v, 10);
		DrawString(img, width, height, x, y, value, borderWidth, v, borderWidth, size, times);	
	}

	template <typename ImgT, typename TValue>
	void DrawString(ImgT * img, int width, int height, int x, int y, TValue value, TValue borderValue, char* text, int borderWidth = 1, CharSizeType size = CharSizeType::_5x7, int times = 1)
	{
		for(int dy = -borderWidth; dy <= borderWidth; ++dy)
			for(int dx = -borderWidth; dx <= borderWidth; ++dx)
				if (dx != 0 && dy != 0)
					DrawString(img, width, height, x + dx, y + dy, borderValue, text, size, times);
		DrawString(img, width, height, x, y, value, text, size, times);
	}

	string toLower(const string & str);
}