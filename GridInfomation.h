#pragma once
#include <cmath>
#include "Vector2D.h"
#include "Array2D.h"

using namespace std;
class GridInfo
{
public:
	union
	{
		struct { int iRes, jRes; };
		int res[2];
	};

	union
	{
		struct { int iStart, jStart, iEnd, jEnd; };
		struct { int ijStart[2], ijEnd[2]; };
	};

	union
	{
		struct { double xMin, yMin, xMax, yMax; };
		struct { double xyMin[2], xyMax[2]; };
	};

	union
	{
		struct { double dx, dy; };
		double dxdy[2];
	};

	union // 2dx, 2dy
	{
		struct { double twodx, twody; };
		double twodxdy[2];
	};

	union // dx^2, dy^2
	{
		struct { double dx2, dy2; };
		double dxdy2[2];
	};

	union // 1/dx, 1/dy
	{
		struct { double oneOverdx, oneOverdy; };
		double oneOverdxdy[2];
	};

	union // 1/2dx, 1/2dy
	{
		struct { double oneOver2dx, oneOver2dy; };
		double oneOverdxdy[2];
	};

	union // 1/dx^2, 1/dy^2
	{
		struct { double oneOverdx2, oneOverdy2; };
		double oneOverdxdy2[2];
	};

	double X0, X1, Y0, Y1, Z0, Z1;
	double deltaX, deltaY, deltaZ;
	int numX, numY, numZ;
	int dimension;

	int numMatX, numMatY, numMatZ;

	double* x;
	double* y;
	double* z;

	GridInfo();
	GridInfo(double x0, double x1, int num1);
	GridInfo(double x0, double x1, int num1, double y0, double y1, int num2);
	GridInfo(double x0, double x1, int num1, double y0, double y1, int num2, double z0, double z1, int num3);
	GridInfo(GridInfo& inputGridInfo);
	
	~GridInfo();

	void gridInitialization(double x0, double x1, int num1);
	void gridInitialization(double x0, double x1, int num1, double y0, double y1, int num2);
	void gridInitialization(double x0, double x1, int num1, double y0, double y1, int num2, double z0, double z1, int num3);

	int index(int i, int j);
	int index(int i, int j, int k);

	GridInfo& operator=(const GridInfo& originGrid)
	{
		if (originGrid.dimension==1)
		{
			gridInitialization(originGrid.X0, originGrid.X1, originGrid.numX);
		}
		if (originGrid.dimension==2)
		{
			gridInitialization(originGrid.X0, originGrid.X1, originGrid.numX, originGrid.Y0, originGrid.Y1, originGrid.numY);
		}
		if (originGrid.dimension==3)
		{
			gridInitialization(originGrid.X0, originGrid.X1, originGrid.numX, originGrid.Y0, originGrid.Y1, originGrid.numY, originGrid.Z0, originGrid.Z1, originGrid.numZ);
		}
		return *this;
	}
private:

};

inline GridInfo::GridInfo()
{
	X0 = 0, X1 = 0, Y0 = 0, Y1 = 0, Z0 = 0, Z1 = 0;
	deltaX = 0, deltaY = 0, deltaZ = 0;
	numX = 0, numY = 0, numZ = 0;
	dimension = 0;

	numMatX = 0, numMatY = 0, numMatZ = 0;

	x = nullptr;
	y = nullptr;
	z = nullptr;
}

GridInfo::GridInfo(double x0, double x1, int num1)
{
	gridInitialization(x0, x1, num1);
}
GridInfo::GridInfo(double x0, double x1, int num1, double y0, double y1, int num2)
{
	gridInitialization(x0, x1, num1, y0, y1, num2);
}
GridInfo::GridInfo(double x0, double x1, int num1, double y0, double y1, int num2, double z0, double z1, int num3)
{
	gridInitialization(x0, x1, num1, y0, y1, num2, z0, z1, num3);
}

GridInfo::GridInfo(GridInfo & inputGridInfo)
{
	if (inputGridInfo.dimension == 1)
	{
		gridInitialization(inputGridInfo.X0, inputGridInfo.X1, inputGridInfo.numX);
	}
	if (inputGridInfo.dimension == 2)
	{
		gridInitialization(inputGridInfo.X0, inputGridInfo.X1, inputGridInfo.numX, inputGridInfo.Y0, inputGridInfo.Y1, inputGridInfo.numY);
	}
	if (inputGridInfo.dimension == 3)
	{
		gridInitialization(inputGridInfo.X0, inputGridInfo.X1, inputGridInfo.numX, inputGridInfo.Y0, inputGridInfo.Y1, inputGridInfo.numY, inputGridInfo.Z0, inputGridInfo.Z1, inputGridInfo.numZ);
	}
}

GridInfo::~GridInfo()
{
	delete[] x, y, z;
}

void GridInfo::gridInitialization(double x0, double x1, int num1)
{
	X0 = x0;
	X1 = x1;
	numX = num1;
	numMatX = numX - 2;
	deltaX = (X1 - X0) / double(numX - 1);
	x = new double[numX];
	dimension = 1;

	for (int i = 0; i < numX; i++)
	{
		x[i] = X0 + deltaX*double(i);
		//cout << x[i] << "\n";
	}
	y = NULL;
	z = NULL;
}

void GridInfo::gridInitialization(double x0, double x1, int num1, double y0, double y1, int num2)
{
	X0 = x0;
	X1 = x1;
	numX = num1;
	numMatX = numX - 2;
	Y0 = y0;
	Y1 = y1;
	numY = num2;
	numMatY = numY - 2;
	deltaX = (X1 - X0) / double(numX - 1);
	deltaY = (Y1 - Y0) / double(numY - 1);
	x = new double[numX];
	y = new double[numY];
	dimension = 2;

	for (int i = 0; i < numX; i++)
	{
		x[i] = X0 + deltaX*double(i);
	}
	for (int i = 0; i < numY; i++)
	{
		y[i] = Y0 + deltaY*double(i);
	}
	z = NULL;
}
void GridInfo::gridInitialization(double x0, double x1, int num1, double y0, double y1, int num2, double z0, double z1, int num3)
{

	X0 = x0;
	X1 = x1;
	numX = num1;
	numMatX = numX - 2;
	Y0 = y0;
	Y1 = y1;
	numY = num2;
	numMatY = numY - 2;
	Z0 = z0;
	Z1 = z1;
	numZ = num3;
	numMatZ = numZ - 2;
	deltaX = (X1 - X0) / double(numX - 1);
	deltaY = (Y1 - Y0) / double(numY - 1);
	deltaZ = (Z1 - Z0) / double(numZ - 1);
	x = new double[numX];
	y = new double[numY];
	z = new double[numZ];
	dimension = 3;

	for (int i = 0; i < numX; i++)
	{
		x[i] = X0 + deltaX*double(i);
	}
	for (int i = 0; i < numY; i++)
	{
		y[i] = Y0 + deltaY*double(i);
	}
	for (int i = 0; i < numZ; i++)
	{
		z[i] = Z0 + deltaZ*double(i);
	}
}

inline int GridInfo::index(int i, int j)
{
	return i + j*numX;
}

inline int GridInfo::index(int i, int j, int k)
{
	return i + j*numX + k*numX*numY;
}

