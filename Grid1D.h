#pragma once

//#ifndef Grid1D_H
//#define Grid1D_H
#include "VectorND.h"


class Grid1D
{
public:
	int iRes;
	int iStart, iEnd;

	double xMin, xMax;

	double xLength;

	double dx;
	double twodx;

	// dx^2
	double dx2;

	// 1/dx
	double oneOverdx;

	// 1/2dx
	double oneOver2dx;

	// 1/dx^2 
	double oneOverdx2;


	Grid1D();
	~Grid1D();

	Grid1D(const double& ipXMin, const double& ipXmax, const int& ipiRes);
	Grid1D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes);
	Grid1D(const Grid1D& ipGrid);

	void initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes);

	inline void operator=(const Grid1D& ipGrid);

	inline double operator ()(const int& i)const;

	//inline Vector2D<double> operator ()(const Vector2D<int> ipVector)const;

	double point(const int& i);
	double cellCenter(const int& i);
	int cellIndex(const double& x);

	inline void Variable();
	inline void Variable(const char * varName1);
private:

};


//#endif // !Grid1D





Grid1D::Grid1D()
{
}

Grid1D::~Grid1D()
{
}

inline Grid1D::Grid1D(const double & ipXMin, const double & ipXmax, const int & ipiRes)
{
	initialize(ipXMin, ipXmax, 0, ipiRes);
}

inline Grid1D::Grid1D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes)
{
	initialize(ipXMin, ipXmax, ipiStart, ipiRes);
}

inline Grid1D::Grid1D(const Grid1D & ipGrid)
{
	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes);
}

inline void Grid1D::initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes)
{
	iRes = ipiRes;
	iStart = ipiStart;
	iEnd = iStart + iRes - 1;
	xMin = ipXMin;
	xMax = ipXmax;
	xLength = xMax - xMin;
	dx = xLength / (double)(iRes - 1);
	twodx = 2.0 * dx;
	dx2 = dx*dx;
	oneOverdx = 1.0 / dx;
	oneOver2dx = 1.0 / twodx;
	oneOverdx2 = oneOverdx*oneOverdx;
}

inline void Grid1D::operator=(const Grid1D & ipGrid)
{
	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes);
}

inline double Grid1D::operator()(const int & i) const
{
	//assert(i >= iStart && i <= iEnd);
	//assert(j >= jStart && j <= jEnd);
	return (xMin + double(i - iStart)*dx);
}

//inline Vector2D<double> Grid1D::operator()(const Vector2D<int> ipVector) const
//{
//	//assert(ipVector.i >= iStart && ipVector.i <= iEnd);
//	//assert(ipVector.j >= jStart && ipVector.j <= jEnd);
//
//	return Vector2D<double>(xMin + double(ipVector.i - iStart)*dx, yMin + double(ipVector.j - jStart)*dy);
//}

inline double Grid1D::point(const int & i)
{
	//assert(i >= iStart && i <= iEnd);
	//assert(j >= jStart && j <= jEnd);

	return (xMin + double(i - iStart)*dx);
}

inline double Grid1D::cellCenter(const int & i)
{
	//assert(i >= iStart && i <= iEnd-1);
	//assert(j >= jStart && j <= jEnd-1);

	return (xMin + (double(i - iStart) + 0.5)*dx);
}

inline int Grid1D::cellIndex(const double & x)
{
	return (floor((x - xMin) + oneOverdx));
}

inline void Grid1D::Variable()
{

	string str = "X=" + to_string(xMin) + ":" + to_string(dx) + ":" + to_string(xMax) + ";";
	const char* cmd = str.c_str();
	MATLAB.Command(cmd);
}

inline void Grid1D::Variable(const char * varName1)
{
	string str = string(varName1) + "=" + to_string(xMin) + ":" + to_string(dx) + ":" + to_string(xMax) + ";";
	const char* cmd = str.c_str();
	MATLAB.Command(cmd);
}

inline std::ostream& operator<<(std::ostream& output, const Grid1D& grid)
{
	output << "GRID_STRUCTURE_2D" << endl;
	output << "- Resolution = " << grid.iRes << endl;
	output << "- Index range =" << grid.iStart << " to " << grid.iEnd << endl;
	output << "- Range = " << grid.xMin << " to " << grid.xMax << endl;
	output << "- dx = " << grid.dx << endl;

	return output;
}

