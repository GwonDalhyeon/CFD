#pragma once
#include <vector>

#include "CommonDef.h"
#include "Vector2D.h"
#include "VectorND.h"
#include "Field2D.h"
#include "Grid2D.h"

class MyBitmap {
public:
	MyBitmap() : pen_(nullptr) {}
	~MyBitmap();


	bool Create(int w, int h);

	void SetPenColor(DWORD clr);

	bool SaveBitmap(const char* path);

	HDC hdc() { return hdc_; }
	int width() { return width_; }
	int height() { return height_; }

private:
	HBITMAP bmp_;
	HDC hdc_;
	HPEN pen_;
	int width_, height_;
};

MyBitmap::~MyBitmap()
{
	DeleteObject(pen_);
	DeleteDC(hdc_);
	DeleteObject(bmp_);

}

inline bool MyBitmap::Create(int w, int h)
{
	BITMAPINFO	bi;
	ZeroMemory(&bi, sizeof(bi));

	bi.bmiHeader.biSize = sizeof(bi.bmiHeader);
	bi.bmiHeader.biBitCount = sizeof(DWORD) * 8;
	bi.bmiHeader.biCompression = BI_RGB;
	bi.bmiHeader.biPlanes = 1;
	bi.bmiHeader.biWidth = w;
	bi.bmiHeader.biHeight = -h;

	void *bits_ptr = nullptr;
	HDC dc = GetDC(GetConsoleWindow());
	bmp_ = CreateDIBSection(dc, &bi, DIB_RGB_COLORS, &bits_ptr, nullptr, 0);
	if (!bmp_) return false;

	hdc_ = CreateCompatibleDC(dc);
	SelectObject(hdc_, bmp_);
	ReleaseDC(GetConsoleWindow(), dc);

	width_ = w;
	height_ = h;

	return true;
}

inline void MyBitmap::SetPenColor(DWORD clr)
{
	if (pen_) DeleteObject(pen_);
	pen_ = CreatePen(PS_SOLID, 1, clr);
	SelectObject(hdc_, pen_);
}

inline bool MyBitmap::SaveBitmap(const char * path)
{
	HANDLE file = CreateFile(path, GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
	if (file == INVALID_HANDLE_VALUE)
	{
		return false;
	}

	BITMAPFILEHEADER fileheader;
	BITMAPINFO infoheader;
	BITMAP bitmap;
	GetObject(bmp_, sizeof(bitmap), &bitmap);

	DWORD* dwp_bits = new DWORD[bitmap.bmWidth * bitmap.bmHeight];
	ZeroMemory(dwp_bits, bitmap.bmWidth * bitmap.bmHeight * sizeof(DWORD));
	ZeroMemory(&infoheader, sizeof(BITMAPINFO));
	ZeroMemory(&fileheader, sizeof(BITMAPFILEHEADER));

	infoheader.bmiHeader.biBitCount = sizeof(DWORD) * 8;
	infoheader.bmiHeader.biCompression = BI_RGB;
	infoheader.bmiHeader.biPlanes = 1;
	infoheader.bmiHeader.biSize = sizeof(infoheader.bmiHeader);
	infoheader.bmiHeader.biHeight = bitmap.bmHeight;
	infoheader.bmiHeader.biWidth = bitmap.bmWidth;
	infoheader.bmiHeader.biSizeImage = bitmap.bmWidth * bitmap.bmHeight * sizeof(DWORD);

	fileheader.bfType = 0x4D42;
	fileheader.bfOffBits = sizeof(infoheader.bmiHeader) + sizeof(BITMAPFILEHEADER);
	fileheader.bfSize = fileheader.bfOffBits + infoheader.bmiHeader.biSizeImage;

	GetDIBits(hdc_, bmp_, 0, height_, (LPVOID)dwp_bits, &infoheader, DIB_RGB_COLORS);

	DWORD wb;
	WriteFile(file, &fileheader, sizeof(BITMAPFILEHEADER), &wb, nullptr);
	WriteFile(file, &infoheader.bmiHeader, sizeof(infoheader.bmiHeader), &wb, nullptr);
	WriteFile(file, dwp_bits, bitmap.bmWidth * bitmap.bmHeight * 4, &wb, nullptr);
	CloseHandle(file);

	delete[] dwp_bits;
	return true;
}
