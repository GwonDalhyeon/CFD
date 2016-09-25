#pragma once


//#ifndef LevelSet2D_H
//#define LevelSet2D_H
#include "CommonDef.h"
#include "Vector2D.h"
#include "VectorND.h"
#include "Grid2D.h"
#include "Field2D.h"



class LevelSet2D
{
public:
	Grid2D grid;

	FD phi;
	FV normal;
	FV unitNormal;
	FV tangential;
	FD meanCurvature;


	// Local Level Set Variables
	Array2D<int> tube;
	VectorND<VI> tubeIndex;
	int numTube;
	double gamma1;
	double gamma2;
	double gamma3;

	LevelSet2D();
	~LevelSet2D();

	LevelSet2D(const Grid2D& ipGrid);
	LevelSet2D(const Grid2D& ipGrid, const double& ipGamma);
	LevelSet2D(const Grid2D& ipGrid, const double& ipGamma1, const double& ipGamma2, const double& ipGamma3);
	LevelSet2D(const LevelSet2D& ipLevelSet, const double& ipGamma);
	LevelSet2D(const LevelSet2D& ipLevelSet, const double& ipGamma1, const double& ipGamma2, const double& ipGamma3);
	LevelSet2D(const FD& ipField);

	const int index(const int& i, const int& j) const;

	const int index(const VI& ipVector) const;

	inline double& operator ()(const int& i, const int& j) const;

	inline double& operator ()(const VI& ipVector) const;

	inline double operator ()(const double& x, const double& y) const;

	inline double operator ()(const VT& ipVector) const;

	inline void operator = (const LevelSet2D& ipLevelSet);

	inline void ComputeNormal();
	inline void LComputeNormal();
	inline VT ComputeNormal(const int& i, const int& j);
	inline VT ComputeNormal(const VI ipVector);

	inline void ComputeUnitNormal();
	inline VT ComputeUnitNormal(const int& i, const int& j);
	inline VT ComputeUnitNormal(const VI ipVector);

	inline void ComputeMeanCurvature();
	inline double ComputeMeanCurvature(const int& i, const int& j);
	inline double ComputeMeanCurvature(const VI & ipVector);

	inline VT gradient(const int& i, const int& j);

	inline double interpolation(const double& x, const double& y);
	inline double interpolation(const VT& ipVector);
	inline double interpolation(const int& i, const int& j);
	//inline double interpolation(const VI& ipVector);

	inline void FillGhostCell();

	// Derivative
	inline double dxxPhi(const int& i, const int& j);
	inline double dxPhi(const int& i, const int& j);
	inline double dxPlusPhi(const int& i, const int& j);
	inline double dxMinusPhi(const int& i, const int& j);

	inline double dyyPhi(const int& i, const int& j);
	inline double dyPhi(const int& i, const int& j);
	inline double dyPlusPhi(const int& i, const int& j);
	inline double dyMinusPhi(const int& i, const int& j);

	inline double dxyPhi(const int& i, const int& j);


	// Local Level Set Functions
	inline double Cutoff(const double & constant);
	inline double Cutoff(const int& i, const int& j);
	inline double Cutoff23(const double & constant);
	inline double Cutoff23(const int& i, const int& j);

	inline void InitialTube();
	inline void UpdateInterface();
	inline void UpdateInterface(const double& ipGamma);
	inline void UpdateInterface(const double& ipGamma1, const double& ipGamma2, const double& ipGamma3);

	inline void UpdateLLS();

	inline void TubeIndex(const int&k, int&i, int& j);
private:

};


//#endif // !LevelSet2D



LevelSet2D::LevelSet2D()
{
}

LevelSet2D::~LevelSet2D()
{
}


inline LevelSet2D::LevelSet2D(const Grid2D & ipGrid)
{
	grid = ipGrid;
	phi = FD(ipGrid);
	normal = FV(ipGrid);
	unitNormal = FV(ipGrid);
	tangential = FV(ipGrid);
	meanCurvature = FD(ipGrid);

	tube = Array2D<int>(ipGrid);

	tubeIndex = VectorND<VI>(1, grid.iRes*grid.jRes);

	gamma1 = 3.0*max(ipGrid.dx, ipGrid.dy);
	gamma2 = 2.0*gamma1;
	gamma3 = 3.0*gamma1;
}

inline LevelSet2D::LevelSet2D(const Grid2D & ipGrid, const double & ipGamma)
{
	grid = ipGrid;
	phi = FD(ipGrid);
	normal = FV(ipGrid);
	unitNormal = FV(ipGrid);
	tangential = FV(ipGrid);
	meanCurvature = FD(ipGrid);

	tube = Array2D<int>(ipGrid);

	tubeIndex = VectorND<VI>(1, grid.iRes*grid.jRes);

	gamma1 = ipGamma;
	gamma2 = 2.0*gamma1;
	gamma3 = 3.0*gamma1;
}

inline LevelSet2D::LevelSet2D(const Grid2D & ipGrid, const double & ipGamma1, const double & ipGamma2, const double & ipGamma3)
{
	grid = ipGrid;
	phi = FD(ipGrid);
	normal = FV(ipGrid);
	unitNormal = FV(ipGrid);
	tangential = FV(ipGrid);
	meanCurvature = FD(ipGrid);

	tube = Array2D<int>(ipGrid);
	tubeIndex = VectorND<VI>(1, grid.iRes*grid.jRes);
	gamma1 = ipGamma1;
	gamma2 = 2.0*gamma1;
	gamma3 = 3.0*gamma1;
}

inline LevelSet2D::LevelSet2D(const LevelSet2D & ipLevelSet, const double & ipGamma)
{
	grid = ipLevelSet.grid;
	phi = FD(grid);
	normal = FV(grid);
	unitNormal = FV(grid);
	tangential = FV(grid);
	meanCurvature = FD(grid);

	tube = Array2D<int>(grid);
	tubeIndex = VectorND<VI>(1, grid.iRes*grid.jRes);
	gamma1 = ipGamma;
	gamma2 = 2.0*gamma1;
	gamma3 = 3.0*gamma1;
}

inline LevelSet2D::LevelSet2D(const LevelSet2D & ipLevelSet, const double & ipGamma1, const double & ipGamma2, const double & ipGamma3)
{
	grid = ipLevelSet.grid;
	phi = FD(grid);
	normal = FV(grid);
	unitNormal = FV(grid);
	tangential = FV(grid);
	meanCurvature = FD(grid);

	tube = Array2D<int>(grid);
	tubeIndex = VectorND<VI>(1, grid.iRes*grid.jRes);
	gamma1 = ipGamma1;
	gamma2 = ipGamma2;
	gamma3 = ipGamma3;
}

inline LevelSet2D::LevelSet2D(const FD& ipField)
{
}

const int LevelSet2D::index(const int & i, const int & j) const
{
	assert(i >= phi.iStart && i <= phi.iEnd);
	assert(j >= phi.jStart && j <= phi.jEnd);
	return phi.index(i, j);
}

const int LevelSet2D::index(const VI& ipVector) const
{
	assert(ipVector[0] >= phi.iStart && ipVector[0] <= phi.iEnd);
	assert(ipVector[1] >= phi.jStart && ipVector[1] <= phi.jEnd);
	return phi.index(ipVector);
}

inline double & LevelSet2D::operator()(const int & i, const int & j) const
{
	assert(i >= phi.iStart && i <= phi.iEnd);
	assert(j >= phi.jStart && j <= phi.jEnd);
	return phi(i, j);
}


inline double & LevelSet2D::operator()(const VI& ipVector) const
{
	assert(ipVector[0] >= phi.iStart && ipVector[0] <= phi.iEnd);
	assert(ipVector[1] >= phi.jStart && ipVector[1] <= phi.jEnd);
	return phi(ipVector);
}



inline double LevelSet2D::operator()(const double & x, const double & y) const
{
	assert(x >= grid.xMin && x <= grid.xMax);
	assert(y >= grid.yMin && y <= grid.yMax);

	return phi(x, y);
}

inline double LevelSet2D::operator()(const VT& ipVector) const
{
	assert(ipVector[0] >= grid.xMin && ipVector[0] <= grid.xMax);
	assert(ipVector[1] >= grid.yMin && ipVector[1] <= grid.yMax);
	return phi(ipVector.x, ipVector.y);
}

inline void LevelSet2D::operator=(const LevelSet2D & ipLevelSet)
{
	grid = ipLevelSet.grid;
	phi = ipLevelSet.phi;
	normal = ipLevelSet.normal;
	unitNormal = ipLevelSet.unitNormal;
	tangential = ipLevelSet.tangential;
	meanCurvature = ipLevelSet.meanCurvature;

	tube = ipLevelSet.tube;
	tubeIndex = ipLevelSet.tubeIndex;
	numTube = ipLevelSet.numTube;
	gamma1 = ipLevelSet.gamma1;
	gamma2 = ipLevelSet.gamma2;
	gamma3 = ipLevelSet.gamma3;
}

inline void LevelSet2D::ComputeNormal()
{
#pragma omp parallel for
	for (int j = phi.jStart; j <= phi.jEnd; j++)
	{
		for (int i = phi.iStart; i <= phi.iEnd; i++)
		{
			normal.dataArray(i, j) = ComputeNormal(i, j);
		}
	}
}

inline void LevelSet2D::LComputeNormal()
{
	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= numTube; k++)
	{
		i = tubeIndex(k).i;
		j = tubeIndex(k).j;
		normal.dataArray(i, j) = ComputeNormal(i, j);
	}
}

inline VT LevelSet2D::ComputeNormal(const int & i, const int & j)
{
	VT normal;

	if (j > phi.jStart && j < phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}
	}
	else if (j == phi.jStart)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else if (j == phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else
	{
		assert(j >= phi.jStart && j <= phi.jEnd);
	}

	normal /= sqrt(normal.x * normal.x + normal.y * normal.y + DBL_EPSILON);
	return normal;
}

inline VT LevelSet2D::ComputeNormal(const VI ipVector)
{
	return ComputeNormal(ipVector[0], ipVector[1]);
}



inline void LevelSet2D::ComputeUnitNormal()
{
#pragma omp parallel for
	for (int j = phi.jStart; j <= phi.jEnd; j++)
	{
		for (int i = phi.iStart; i <= phi.iEnd; i++)
		{
			unitNormal.dataArray(i, j) = ComputeUnitNormal(i, j);
		}
	}
}

inline VT LevelSet2D::ComputeUnitNormal(const int & i, const int & j)
{
	VT normal;

	if (j > phi.jStart && j < phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}
	}
	else if (j == phi.jStart)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else if (j == phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else
	{
		assert(j >= phi.jStart && j <= phi.jEnd);
	}


	return normal;
}

inline VT LevelSet2D::ComputeUnitNormal(const VI ipVector)
{
	return ComputeUnitNormal(ipVector[0], ipVector[1]);
}

inline void LevelSet2D::ComputeMeanCurvature()
{
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.iStart; j <= grid.jEnd; j++)
		{
			meanCurvature(i, j) = ComputeMeanCurvature(i, j);
		}
	}
}

inline double LevelSet2D::ComputeMeanCurvature(const int & i, const int & j)
{
	return -(dxxPhi(i, j)*dyPhi(i, j)*dyPhi(i, j) - 2.0*dxyPhi(i, j)*dxPhi(i, j)*dyPhi(i, j) + dyyPhi(i, j)*dxPhi(i, j)*dxPhi(i, j)) / pow(dxPhi(i, j)*dxPhi(i, j) + dyPhi(i, j)*dyPhi(i, j) + DBL_EPSILON, 3.0 / 2.0);
}

inline double LevelSet2D::ComputeMeanCurvature(const VI & ipVector)
{
	return ComputeMeanCurvature(ipVector.i, ipVector.j);
}

inline VT LevelSet2D::gradient(const int & i, const int & j)
{
	return VT(dxPhi(i, j), dyPhi(i, j));
}

inline double LevelSet2D::interpolation(const double & x, const double & y)
{
	if (grid.xMin <= x && grid.xMax >= x &&grid.yMin <= y && grid.yMax >= y)
	{
		VT xy(x, y);
		VI cell = phi.containedCell(x, y);

		double distance00, distance10, distance01, distance11;
		distance00 = (grid.point(cell.i, cell.j) - xy).magnitude();
		distance10 = (grid.point(cell.i + 1, cell.j) - xy).magnitude();
		distance01 = (grid.point(cell.i, cell.j + 1) - xy).magnitude();
		distance11 = (grid.point(cell.i + 1, cell.j + 1) - xy).magnitude();
		if (distance00 < grid.dx / 2)
		{
			return phi(cell);
		}
		if (distance10 < grid.dx / 2)
		{
			return phi(cell.i + 1, cell.j);
		}
		if (distance01 < grid.dx / 2)
		{
			return phi(cell.i, cell.j + 1);
		}
		if (distance11 < grid.dx / 2)
		{
			return phi(cell.i + 1, cell.j + 1);
		}

		return ((phi(cell)*distance00 + phi(cell.i + 1, cell.j)*distance10 + phi(cell.i, cell.j + 1)*distance01 + phi(cell.i + 1, cell.j + 1)*distance11) / (distance00 + distance01 + distance10 + distance11));
	}
	else if (grid.xMin <= x && grid.xMax >= x &&grid.yMin > y)
	{
		double value1 = interpolation(x, grid.yMin);
		double value2 = interpolation(x, grid.yMin + grid.dy);

		return (value1 - value2)*grid.oneOverdy*(grid.yMin - y) + value1;
	}
	else if (grid.xMin <= x && grid.xMax >= x &&grid.yMax < y)
	{
		double value1 = interpolation(x, grid.yMax);
		double value2 = interpolation(x, grid.yMax - grid.dy);

		return (value1 - value2)*grid.oneOverdy*(y - grid.yMax) + value1;
	}
	else if (grid.xMin > x  &&grid.yMin <= y && grid.yMax >= y)
	{
		double value1 = interpolation(grid.xMin, y);
		double value2 = interpolation(grid.xMin + grid.dx, y);

		return (value1 - value2)*grid.oneOverdx*(grid.xMin - x) + value1;
	}
	else if (grid.xMax < x &&grid.yMin <= y && grid.yMax >= y)
	{
		double value1 = interpolation(grid.xMax, y);
		double value2 = interpolation(grid.xMax - grid.dx, y);

		return (value1 - value2)*grid.oneOverdx*(x - grid.xMax) + value1;
	}
	else
	{
		double tempX;
		double tempY;
		if (x < grid.xMin)
		{
			tempX = grid.xMin;
		}
		else
		{
			tempX = grid.xMax;
		}

		if (y < grid.yMin)
		{
			tempY = grid.yMin;
		}
		else
		{
			tempY = grid.yMax;
		}

		double value1 = interpolation(tempX, y);
		double value2 = interpolation(x, tempY);

		return (value1 + value2) / 2.0;
	}
}

inline double LevelSet2D::interpolation(const VT& ipVector)
{
	return interpolation(ipVector.x, ipVector.y);
}

inline double LevelSet2D::interpolation(const int & i, const int & j)
{
	return interpolation(grid.point(i, j));
}

inline void LevelSet2D::FillGhostCell()
{
	phi.FillGhostCell();
	normal.FillGhostCell();
	unitNormal.FillGhostCell();
	tangential.FillGhostCell();
	meanCurvature.FillGhostCell();
}

inline double LevelSet2D::dxxPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (phi(i + 1, j) - 2.0 * phi(i, j) + phi(i - 1, j))*grid.oneOverdx2;
	}
	else if (i == grid.iStart)
	{
		double tempPhi = interpolation(i - 1, j);
		return (phi(i + 1, j) - 2.0 * phi(i, j) + tempPhi)*grid.oneOverdx2;
		//return (phi(i, j) - 2 * phi(i + 1, j) + phi(i + 2, j))*grid.oneOver2dx;
	}
	else
	{
		double tempPhi = interpolation(i + 1, j);
		return (tempPhi - 2.0 * phi(i, j) + phi(i - 1, j))*grid.oneOverdx2;
		//return (phi(i - 2, j) - 2 * phi(i - 1, j) + phi(i, j))*grid.oneOver2dx;
	}
}

inline double LevelSet2D::dxPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (phi(i + 1, j) - phi(i - 1, j))*grid.oneOver2dx;
	}
	else if (i == grid.iStart)
	{
		double tempPhi = interpolation(i - 1, j);
		return (phi(i + 1, j) - tempPhi)*grid.oneOver2dx;
		//return dxPlusPhi(i, j);
	}
	else
	{
		double tempPhi = interpolation(i + 1, j);
		return (tempPhi - phi(i - 1, j))*grid.oneOver2dx;
		//return dxMinusPhi(i, j);
	}
}

inline double LevelSet2D::dxPlusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i < grid.iEnd)
	{
		return (phi(i + 1, j) - phi(i, j))*grid.oneOverdx;
	}
	else
	{
		double tempPhi = interpolation(i + 1, j);
		return (tempPhi - phi(i, j))*grid.oneOverdx;
		//return (phi(i, j) - phi(i - 1, j))*grid.oneOverdx;
	}
}

inline double LevelSet2D::dxMinusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart)
	{
		return (phi(i, j) - phi(i - 1, j))*grid.oneOverdx;
	}
	else
	{
		double tempPhi = interpolation(i - 1, j);
		return (phi(i, j) - tempPhi)*grid.oneOverdx;
		//return (phi(i + 1, j) - phi(i, j))*grid.oneOverdx;
	}
}

inline double LevelSet2D::dyyPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (phi(i, j + 1) - 2 * phi(i, j) + phi(i, j - 1))*grid.oneOverdy2;
	}
	else if (j == grid.jStart)
	{
		return (phi(i, j) - 2 * phi(i, j + 1) + phi(i, j + 2))*grid.oneOverdy2;
	}
	else
	{
		return (phi(i, j - 2) - 2 * phi(i, j - 1) + phi(i, j))*grid.oneOverdy2;
	}
}

inline double LevelSet2D::dyPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (phi(i, j + 1) - phi(i, j - 1))*grid.oneOver2dy;
	}
	else if (j == grid.jStart)
	{
		return dyPlusPhi(i, j);
	}
	else
	{
		return dyMinusPhi(i, j);
	}
}

inline double LevelSet2D::dyPlusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j < grid.jEnd)
	{
		return (phi(i, j + 1) - phi(i, j))*grid.oneOverdy;
	}
	else
	{
		return (phi(i, j) - phi(i, j - 1))*grid.oneOverdy;
	}
}

inline double LevelSet2D::dyMinusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart)
	{
		return (phi(i, j) - phi(i, j - 1))*grid.oneOverdy;
	}
	else
	{
		return (phi(i, j + 1) - phi(i, j))*grid.oneOverdy;
	}
}

inline double LevelSet2D::dxyPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (dxPhi(i, j + 1) - dxPhi(i, j - 1))*grid.oneOver2dy;
	}
	else if (j == grid.iStart)
	{
		return (dxPhi(i, j + 1) - dxPhi(i, j))*grid.oneOverdy;
	}
	else
	{
		return (dxPhi(i, j) - dxPhi(i, j - 1))*grid.oneOverdy;
	}
}

inline double LevelSet2D::Cutoff(const double & constant)
{
	double absConst = abs(constant);
	if (gamma2 < absConst)
	{
		return 0.0;
	}
	else if (gamma1 < absConst && absConst <= gamma2)
	{
		return (absConst - gamma2)*(absConst - gamma2)*(2.0*absConst + gamma2 - 3.0*gamma1) / pow(gamma2 - gamma1, 3.0);
	}
	else
	{
		return 1.0;
	}
}

inline double LevelSet2D::Cutoff(const int & i, const int & j)
{
	return Cutoff(phi(i, j));
}

inline double LevelSet2D::Cutoff23(const double & constant)
{
	double absConst = abs(constant);
	if (gamma3 < absConst)
	{
		return 0.0;
	}
	else if (gamma2 < absConst && absConst <= gamma3)
	{
		return (absConst - gamma3)*(absConst - gamma3)*(2.0*absConst + gamma3 - 3.0*gamma1) / pow(gamma3 - gamma1, 3.0);
	}
	else
	{
		return 1.0;
	}
}

inline double LevelSet2D::Cutoff23(const int & i, const int & j)
{
	return Cutoff23(phi(i, j));
}

inline void LevelSet2D::InitialTube()
{
	numTube = 0;
	double absConst;
	for (int i = phi.grid.iStart; i <= phi.grid.iEnd; i++)
	{
		for (int j = phi.grid.iStart; j <= phi.grid.jEnd; j++)
		{
			absConst = abs(phi(i, j));
			if (absConst < gamma3)
			{
				tube(i, j) = 3;
				numTube += 1;
				tubeIndex(numTube) = VI(i, j);
				if (absConst < gamma2)
				{
					tube(i, j) = 2;
					if (absConst < gamma1)
					{
						tube(i, j) = 1;
						//numTube += 1;
						//tubeIndex(numTube) = VI(i, j);
					}
				}
			}
		}
	}

	double threshold = gamma2;
	for (int i = phi.grid.iStart; i <= phi.grid.iEnd; i++)
	{
		for (int j = phi.grid.iStart; j <= phi.grid.jEnd; j++)
		{
			if (phi(i, j) < -threshold)
			{
				phi(i, j) = -(threshold);
			}
			else if (phi(i, j) > threshold)
			{
				phi(i, j) = threshold;
			}
		}
	}
}

inline void LevelSet2D::UpdateInterface()
{
	UpdateInterface(gamma1, gamma2, gamma3);
}

inline void LevelSet2D::UpdateInterface(const double & ipGamma)
{
	gamma1 = ipGamma;
	gamma2 = 2.0*ipGamma;
	gamma3 = 3.0*ipGamma;
	UpdateInterface(gamma1, gamma2, gamma3);
}

inline void LevelSet2D::UpdateInterface(const double & ipGamma1, const double & ipGamma2, const double & ipGamma3)
{
	gamma1 = ipGamma1;
	gamma2 = ipGamma2;
	gamma3 = ipGamma3;
	int tempNumTube = numTube;
	tube = 0;
	numTube = 0;
	int i, j;
	double absConst;
	for (int k = 1; k <= tempNumTube; k++)
	{
		i = tubeIndex(k).i;
		j = tubeIndex(k).j;
		absConst = abs(phi(i, j));
		//if (absConst < gamma3 && Tube(i, j) == 0)
		{
			//Tube(i, j) = 3;
			//numTube += 1;
			//tubeIndex(numTube) = VI(i, j);
			if (absConst < gamma2 && tube(i, j) == 0)
			{
				tube(i, j) = 2;
				numTube += 1;
				tubeIndex(numTube) = VI(i, j);
				if (absConst <= gamma1)
				{
					tube(i, j) = 1;
				}
			}
		}
	}
	//Tube.Variable("Tube1");
	//MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	//MATLAB.Command("subplot(1,2,1)");
	//MATLAB.Command("surf(X, Y, Tube1);grid on;axis([-1 1 -1 1]);axis equal;");
	tempNumTube = numTube;
	int bdryWidth = 3;
	for (int k = 1; k <= tempNumTube; k++)
	{
		i = tubeIndex(k).i;
		j = tubeIndex(k).j;
		if (tube(i, j) == 2)
		{
			for (int ii = -bdryWidth; ii <= bdryWidth; ii++)
			{
				for (int jj = -bdryWidth; jj <= bdryWidth; jj++)
				{
					if (tube(i + ii, j + jj) == 0)
					{
						tube(i + ii, j + jj) = 3;
						numTube += 1;
						tubeIndex(numTube) = VI(i + ii, j + jj);
					}
				}
			}
		}
	}

	//Tube.Variable("Tube2");
	//MATLAB.Command("subplot(1,2,2)");
	//MATLAB.Command("surf(X, Y, Tube2);grid on;axis([-1 1 -1 1]);axis equal;");

}

inline void LevelSet2D::UpdateLLS()
{
	double threshold = gamma2;
	int i, j;
	int bdryWidth = 2;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= numTube; k++)
	{
		i = tubeIndex(k).i;
		j = tubeIndex(k).j;
		if (tube(i, j) == 3)
		{
			for (int ii = -bdryWidth; ii <= bdryWidth; ii++)
			{
				for (int jj = -bdryWidth; jj <= bdryWidth; jj++)
				{
					if (tube(i + ii, j + jj) == 0 || tube(i + ii, j + jj) == 3)
					{
						if (phi(i + ii, j + jj) < -threshold)
						{
							phi(i + ii, j + jj) = -(threshold);
						}
						else if (phi(i + ii, j + jj) > threshold)
						{
							phi(i + ii, j + jj) = threshold;
						}
					}
				}
			}
		}
	}
}

inline void LevelSet2D::TubeIndex(const int & k, int & i, int & j)
{
	i = tubeIndex(k).i;
	j = tubeIndex(k).j;
}

typedef LevelSet2D LS;

