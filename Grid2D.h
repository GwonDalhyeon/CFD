#pragma once

//#ifndef Grid2D_H
//#define Grid2D_H
#include "Vector2D.h"


class Grid2D
{
public:
	union
	{
		struct { int iRes, jRes; };
		int ijRes[2];
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
		struct { double xLength, yLength; };
		struct { double xyLength[2]; };
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


	Grid2D();
	~Grid2D();

	Grid2D(const double& ipXMin, const double& ipXmax, const int& ipiRes, const double& ipYMin, const double& ipYmax, const int& ipjRes);
	Grid2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);
	Grid2D(const Grid2D& ipGrid);

	void Initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);

	inline void operator=(const Grid2D& ipGrid);

	inline Vector2D<double> operator ()(const int& i, const int& j)const;

	inline Vector2D<double> operator ()(const Vector2D<int> ipVector)const;

	Vector2D<double> point(const int& i, const int& j);
	Vector2D<double> cellCenter(const int& i, const int& j);
	Vector2D<int> cellIndex(const double& x, const double& y);

	inline void Variable();
	inline void Variable(const char * varName1, const char * varName2);
private:

};


//#endif // !Grid2D





Grid2D::Grid2D()
{
}

Grid2D::~Grid2D()
{
}

inline Grid2D::Grid2D(const double & ipXMin, const double & ipXmax, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjRes)
{
	Initialize(ipXMin, ipXmax, 0, ipiRes, ipYMin, ipYmax, 0, ipjRes);
}

inline Grid2D::Grid2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	Initialize(ipXMin, ipXmax, ipiStart, ipiRes, ipYMin, ipYmax, ipjStart, ipjRes);
}

inline Grid2D::Grid2D(const Grid2D & ipGrid)
{
	Initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.jStart, ipGrid.jRes);
}

inline void Grid2D::Initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	iRes = ipiRes;
	jRes = ipjRes;
	iStart = ipiStart;
	jStart = ipjStart;
	iEnd = iStart + iRes - 1;
	jEnd = jStart + jRes - 1;
	xMin = ipXMin;
	yMin = ipYMin;
	xMax = ipXmax;
	yMax = ipYmax;
	xLength = xMax - xMin;
	yLength = yMax - yMin;
	dx = xLength / (double)(iRes - 1);
	dy = yLength / (double)(jRes - 1);
	twodx = 2.0 * dx;
	twody = 2.0*dy;
	dx2 = dx*dx;
	dy2 = dy*dy;
	oneOverdx = 1.0 / dx;
	oneOverdy = 1.0 / dy;
	oneOver2dx = 1.0 / twodx;
	oneOver2dy = 1.0 / twody;
	oneOverdx2 = oneOverdx*oneOverdx;
	oneOverdy2 = oneOverdy*oneOverdy;
}

inline void Grid2D::operator=(const Grid2D & ipGrid)
{
	Initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.jStart, ipGrid.jRes);
}

inline Vector2D<double> Grid2D::operator()(const int & i, const int & j) const
{
	//assert(i >= iStart && i <= iEnd);
	//assert(j >= jStart && j <= jEnd);
	return Vector2D<double>(xMin + double(i - iStart)*dx, yMin + double(j - jStart)*dy);
}

inline Vector2D<double> Grid2D::operator()(const Vector2D<int> ipVector) const
{
	//assert(ipVector.i >= iStart && ipVector.i <= iEnd);
	//assert(ipVector.j >= jStart && ipVector.j <= jEnd);

	return Vector2D<double>(xMin + double(ipVector.i - iStart)*dx, yMin + double(ipVector.j - jStart)*dy);
}

inline Vector2D<double> Grid2D::point(const int & i, const int & j)
{
	//assert(i >= iStart && i <= iEnd);
	//assert(j >= jStart && j <= jEnd);

	return Vector2D<double>(xMin + double(i - iStart)*dx, yMin + double(j - jStart)*dy);
}

inline Vector2D<double> Grid2D::cellCenter(const int & i, const int & j)
{
	//assert(i >= iStart && i <= iEnd-1);
	//assert(j >= jStart && j <= jEnd-1);

	return Vector2D<double>(xMin + (double(i - iStart) + 0.5)*dx, yMin + (double(j - jStart) + 0.5)*dy);
}

inline Vector2D<int> Grid2D::cellIndex(const double & x, const double & y)
{
	return Vector2D<int>(floor((x - xMin) + oneOverdx), floor((y - yMin) + oneOverdy));
}

inline void Grid2D::Variable()
{
	
	string str = "[X,Y]=meshgrid(" +to_string(xMin) +":" +to_string(dx)+":"+to_string(xMax)+"," + to_string(yMin) + ":" + to_string(dy) + ":" + to_string(yMax) + +");";
	const char* cmd = str.c_str();
	MATLAB.Command(cmd);
}

inline void Grid2D::Variable(const char * varName1, const char * varName2)
{
	string str = "["+string(varName1)+","+ string(varName2) +"]=meshgrid(" + to_string(xMin) + ":" + to_string(dx) + ":" + to_string(xMax) + "," + to_string(yMin) + ":" + to_string(dy) + ":" + to_string(yMax) + +");";
	const char* cmd = str.c_str();
	MATLAB.Command(cmd);
}

inline std::ostream& operator<<(std::ostream& output, const Grid2D& grid)
{
	output << "GRID_STRUCTURE_2D" << endl;
	output << "- Resolution = (" << grid.iRes << "," << grid.jRes << ")" << endl;
	output << "- Index range =(" << grid.iStart << "," << grid.jStart << ") to (" << grid.iEnd << "," << grid.jEnd << ") " << endl;
	output << "- Range = (" << grid.xMin << "," << grid.yMin << ") to (" << grid.xMax << "," << grid.yMax << ") " << endl;
	output << "- dx = (" << grid.dx << "," << grid.dy << ")" << endl;

	return output;
}

